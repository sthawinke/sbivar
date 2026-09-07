#' Calculate bivariate Moran's I between a point pattern and a quantitative modality, with variance and p-value
#'
#' The variance calculation requires estimation of the spatial autocorrelation structure of every feature of the second modality, using Matheron's variogram estimator \insertCite{Matheron1963}{sbivar}.
#'
#' @inheritParams sbivarSinglePPP
#' @param variogramModels A character vector, indicating the variogram model passed onto \link[gstat]{vgm}.
#' Currently, only "Exp" and "Lin" are implemented for computational reasons.
#' @param etas Vectors of weight matrix parameters, whose elements are passed onto \link{buildWeightMat}
#' @param wo type of weight parameter, passed onto \link{buildWeightMat}
#' @param cutoff,width Cutoff and width of the variogram estimation, passed onto \link[gstat]{vgm}
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
#' The maximum value of the bivariate Moran's I statistics are different for every gene pair unlike for \link{MoransISingle} as the weight matrix is random. As this presents too much computation, no maximum values are calculated.
#' @note No multithreading is implemented for the variance calculation, as the matrix calculations involved
#' may use inherent multithreading with OpenBLAS.
#' @importFrom spatstat.geom npoints coords split.ppp coords<-
MoransISinglePPP <- function(X, Y, Ey, wo, etas, cutoff, width, verbose,
    variogramModels, returnSEsMoransI, featuresX, featuresY, findVariances = TRUE, ...) {
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
    movedCoords <- moveTwoCoords(as.matrix(coords(X)), Ey)
    Marks <- marks(X, drop = FALSE)
    coords(X) <- movedCoords$Cx
    Ey <- movedCoords$Ey
    if (verbose) {
        message("Calculating bivariate Moran's I statistics ...")
    }
    wParams <- selfName(switch(wo,
        "Gauss" = etas
    ))
    X <- split.ppp(X, factor(Marks$feature))
    mm2 <- m * (m - 1) / 2
    distY <- as.vector(stats::dist(Ey))
    if (findVariances) {
        # Estimate spatial autocorrelation
        if (verbose) {
            message("Fitting variograms for second modality (", k, " features) ...")
        }
        variogramsY <- matheronVariograms(Y[, featuresY, drop = FALSE], Ey,
            width = width, cutoff = cutoff,
            variogramModels = variogramModels, ...
        )
    }
    res <- lapply(featuresX, function(featx) {
        Cx <- coords(X[[featx]])
        n <- nrow(Cx)
        prodFac <- (n - 1) * (m - 1)
        Ws <- vapply(wParams, FUN.VALUE = matrix(0, n, m), function(iter) {
            buildWeightMat(Cx = Cx, Ey = Ey, wo = wo, eta = iter, numNN = iter)
        })
        Ws <- Ws[, , idW <- (colSums(Ws, dims = 2, na.rm = TRUE) > 0), drop = FALSE]
        numWs <- dim(Ws)[3]
        if (!all(idW) && (wo == "Gauss")) {
            etas <- etas[idW]
        }
        Ixys <- t(t(vapply(seq_len(numWs), FUN.VALUE = double(k), function(i) {
            colSums(Ws[, , i] %*% Y[, featuresY, drop = FALSE])
        })) / sqrt(prodFac)) # Normalize for matrix size
        out <- if (findVariances) {
            varIxy <- t(vapply(selfName(featuresY), FUN.VALUE = double(numWs), function(featy) {
                # C++: build Sigma_X and batch-compute t(W[,,i]) Sigma_X W[,,i] for all i,
                # returning their traces
                sigRes <- computeSigXws(evalVariogram(variogramsY[[featy]], distY), Ws, findSigXws = FALSE)
                # Precomputing evalVariogram for all Y's is too much memory, so repeat it at a speed cost
                return(sigRes$traces)
            }))
            for (i in seq_len(numWs)) { # If negative variance, fall back on independence
                if (length(zeroId <- c(which(varIxy[, i] <= 0), which(is.na(varIxy[, i]))))) {
                    varIxy[zeroId, i] <- sum(Ws[, , i]^2) # tr(W^tW)
                }
            }
            list("Ixys" = Ixys, "seIxy" = sqrt(varIxy / prodFac))
        } else {
            list("Ixys" = Ixys)
        }
        printProgress(featx, featuresX, verbose)
        return(out)
    })
    Ixy <- do.call(what = rbind, lapply(res, function(x) x$Ixys))
    colnames(Ixy) <- paste0("Ixy_", wParams)
    if (findVariances) {
        seIxy <- do.call(what = rbind, lapply(res, function(x) x$seIxy))
        colnames(seIxy) <- paste0("SE(Ixy)_", wParams)
        IxyPvals <- makePval(Ixy / seIxy)
        cctPvals <- apply(IxyPvals, 1, CCT) # CCT correction
        out <- cbind(Ixy, if (returnSEsMoransI) seIxy, "pVal" = cctPvals)
    } else {
        out <- Ixy
    }
    rownames(out) <- make.names(apply(expand.grid(featuresY, featuresX)[, 2:1], 1, paste, collapse = "__"))
    return(list(
        "res" = out
    ))
}
