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
#' @importFrom smoppix loadBalanceBplapply
MoransISinglePPP <- function(X, Y, Ey, wo, etas, cutoff, width, verbose,
    variogramModels, returnSEsMoransI, featuresX, featuresY, findVariances = TRUE, ...) {
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
    featsVec <- marks(X, drop = FALSE)$feature
    wParams <- selfName(switch(wo,
        "Gauss" = etas
    ))
    if (any(wParams > 2e-4)) {
        warning("Eta values larger than 2e-3 are not meaningful for point patterns and may lead to false positive findings!")
    }
    if (findVariances) {
        if (verbose) {
            message("Fitting variograms for second modality (", k, " features) ...")
        }
        variogramsY <- matheronVariograms(Y[, featuresY, drop = FALSE], movedCoords$Ey,
            width = width, cutoff = cutoff,
            variogramModels = variogramModels, ...
        )
        # Pack variogram parameters into a k x 3 matrix for C++.
        # distY is not pre-computed here; it is computed inside C++ per
        # iteration and freed on return.
        vgParY <- do.call(rbind, lapply(featuresY, function(fy) {
            vg <- variogramsY[[fy]]
            c(vg[2L, "psill"], vg[2L, "range"], as.numeric(vg[2L, "model"] == "Exp"))
        }))
    }
    # Prepare a vgParY argument safe to pass regardless of findVariances
    vgParY_arg <- if (findVariances) vgParY else matrix(0.0, 0L, 3L)
    if (verbose) {
        message("Calculating bivariate Moran's I statistics and variances ...")
    }
    res <- loadBalanceBplapply(featuresX, function(featx) {
        Cx_i <- movedCoords$Cx[featsVec == featx, , drop = FALSE]
        n <- nrow(Cx_i)
        prodFac <- (n - 1) * (m - 1)
        # Loop over weight parameters; W is built inside C++ and never
        # materialised in the R session.
        nW <- length(wParams)
        IxysList_w <- vector("list", nW)
        varIxyList_w <- if (findVariances) vector("list", nW)
        idW <- logical(nW)
        for (wi in seq_len(nW)) {
            res_w <- computeIxyAndTracePPP_cpp(
                Cx_i, movedCoords$Ey, wParams[[wi]],
                Y[, featuresY, drop = FALSE],
                vgParY_arg, sqrt(prodFac), findVariances
            )
            if (res_w$isZero) next
            idW[wi] <- TRUE
            IxysList_w[[wi]] <- res_w$Ixys
            if (findVariances) {
                traces <- res_w$traces
                traces[traces <= 0 | is.na(traces)] <- res_w$trWtW
                varIxyList_w[[wi]] <- traces
            }
        }
        if (!all(idW) && (wo == "Gauss")) {
            etas <- etas[idW] # local shadow, preserving original behaviour
        }
        IxysList_w <- IxysList_w[idW]
        if (findVariances) varIxyList_w <- varIxyList_w[idW]
        numWs <- sum(idW)
        # Assemble k x numWs matrices from per-W k-vectors
        Ixys <- do.call(cbind, IxysList_w)
        out <- if (findVariances) {
            varIxy <- do.call(cbind, varIxyList_w) # k x numWs
            list("Ixys" = Ixys, "seIxy" = sqrt(varIxy / prodFac))
        } else {
            list("Ixys" = Ixys)
        }
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
