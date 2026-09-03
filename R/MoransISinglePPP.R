#' Calculate bivariate Moran's I between two modality matrix, with variance and p-value
#'
#' The variance calculation requires estimation of the spatial autocorrelation structure of every feature separately, using Matheron's variogram estimator \insertCite{Matheron1963}{sbivar}.
#'
#' @inheritParams sbivarSinglePPP
#' @param variogramModels A character vector, indicating the variogram model passed onto \link[gstat]{vgm}.
#' Currently, only "Exp" and "Lin" are implemented for computational reasons.
#' @param numNNs,etas Vectors of weight matrix parameters, whose elements are passed onto \link{buildWeightMat}
#' @param wo type of weight parameter, passed onto \link{buildWeightMat}
#' @param cutoff,width Cutoff and width of the variogram estimation, passed onto \link[gstat]{vgm}
#' @param findMaxW Is the maximum bivariate Moran's I needed?
#' @param returnSEsMoransI A boolean, are standard errors of Moran's I to be returned?
#' @param findVariances Should variances be calculated? For internal use only
#' @param ... passed onto \link[gstat]{variogram}
#'
#' @returns A dataframe of results sorted by p-value, also containing the estimated Moran's I statistic and its variance.
#' In addition, the maximum value of the Moran's I statistic, and the parameters of the weight matrix
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats dist
#'
#' @details By default, a number of range parameters and corresponding weight matrices are screened for spatial association,
#' and their p-value combined using the Cauchy combination rule by \insertCite{Liu2020}{sbivar}.
#' The maximum value of the bivariate Moran's I statistics are returned conditionally,
#' as it is computation intensive and not always needed.
#' @note No multithreading is implemented for the variance calculation, as the matrix calculations involved
#' may use inherent multithreading with OpenBLAS.
#' @importFrom spatstat.geom npoints coords split.ppp
MoransISinglePPP <- function(
      X, Y, Ey, wo, etas, numNNs, cutoff, width, verbose,
      findMaxW, variogramModels, returnSEsMoransI, featuresX, featuresY, findVariances = TRUE, ...
) {
    n <- npoints(X)
    m <- nrow(Y)
    p <- length(featuresX)
    k <- length(featuresY)
    if (verbose) {
        message("Testing significance of bivariate Moran's I for ", p * k, " feature pairs")
    }
    # Scale outcomes
    Y <- scale(Y)
    # Move coordinates
    movedCoords <- moveTwoCoords(coords(X), Ey)
    coords(X) <- movedCoords$Cx
    Ey <- movedCoords$Ey
    prodFac <- (n - 1) * (m - 1)
    if (verbose) {
        message("Calculating bivariate Moran's I statistics ...")
    }
    wParams <- selfName(switch(wo,
        "Gauss" = etas,
        "nn" = numNNs
    ))
    X <- split.ppp(X, marks(X, drop = FALSE)$features)
    res <- lapply(featuresX, function(featx) {
        Cx <- coords(X[[featx]])
        Ws <- vapply(wParams, FUN.VALUE = matrix(0, n, m), function(iter) {
            buildWeightMat(Cx = Cx, Ey = Ey, wo = wo, eta = iter, numNN = iter)
        })
        Ws <- Ws[, , idW <- (colSums(Ws, dims = 2, na.rm = TRUE) > 0), drop = FALSE]
        numWs <- dim(Ws)[3]
        if (!all(idW) && (wo == "Gauss")) {
            etas <- etas[idW]
        }
        Ixys <- vapply(seq_len(numWs), FUN.VALUE = matrix(0, p, k), function(i) {
            rowSums(crossprod(Ws[, , i] %*% Y[, featuresY, drop = FALSE]))
        }) / sqrt(prodFac) # Normalize for matrix size
        # Reformat to long format
        out <- matrix(c(Ixys), ncol = numWs, dimnames = list(NULL, paste0("Ixy_", wParams)))
        if (findVariances) {
            # Estimate spatial autocorrelation
            if (verbose) {
                message("Fitting variograms for second modality (", k, " features) ...")
            }
            variogramsY <- matheronVariograms(Y[, featuresY, drop = FALSE], Ey,
                width = width, cutoff = cutoff,
                variogramModels = variogramModels, ...
            )
            distX <- as.vector(stats::dist(Cx))
            distY <- as.vector(stats::dist(Ey))
            if (verbose) {
                message("Calculating variances of bivariate Moran's I statistics ...")
            }
            mm2 <- m * (m - 1) / 2
            varIxy <- vapply(selfName(featuresY), FUN.VALUE = double(numWs), function(featy) {
                # C++: build Sigma_X and batch-compute t(W[,,i]) Sigma_X W[,,i] for all i,
                # returning lower-triangle columns (sigXws, mm2 x numWs) and traces
                sigRes <- computeSigXws(evalVariogram(variogramsY[[featy]], distY), Ws)
                # Precomputing evalVariogram for all Y's is too much memory, so repeat it at a speed cost
                return(sigRes$traces)
            })
            varIxy <- aperm(varIxy, perm = 3:1) # Rearrange
            for (i in seq_len(numWs)) { # If negative variance, fall back on independence
                if (length(zeroId <- c(which(varIxy[, , i] <= 0), which(is.na(varIxy[, , i]))))) {
                    varIxy[, , i][zeroId] <- sum(Ws[, , i]^2) # tr(W^tW)
                }
            }
            varIxy <- varIxy / prodFac # Correct for matrix size
        }
        printProgress(featx, featuresX, verbose)
    })
    # P-values
    IxyPvals <- makePval(Ixys / (seIxy <- sqrt(varIxy)))
    # CCT correction
    cctPvals <- apply(IxyPvals, c(1, 2), CCT)
    if (returnSEsMoransI) {
        out <- cbind(out, matrix(c(seIxy),
            ncol = numWs,
            dimnames = list(NULL, paste0("SE(Ixy)_", wParams))
        ))
    }
    out <- cbind(out, "pVal" = c(cctPvals))
    rownames(out) <- makeNames(featuresX, featuresY)
    return(list(
        "res" = res, "maxIxy" = maxIxy
    ))
}
