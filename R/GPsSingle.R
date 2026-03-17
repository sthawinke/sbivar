#' Fit Gaussian processes (GPs) if needed, and perform score tests
#'
#' Fit all univariate GPs on both modalities, and perform all bivariate tests across them
#'
#' @inheritParams sbivarSingle
#' @inheritParams fitGP
#' @param gpParams Parameters of the Gaussian processes, see details
#' @param numLscAlts Number of length scales to be tested for bivariate association
#' @param Quants Most extreme quantiles of the distance distribution to take as length scales
#' @returns A named list of results
#' @details gpParams must be a list of length 2 with names 'X' and 'Y', consisting of matrices
#' with rownames "mean", "nugget", "range" and "sigma", and column names as in X and Y.
#' This argument allows to pass parameters of the Gaussian processes estimated with other software
#' (e.g. with GPU acceleration) to perform the score test.
GPsSingle <- function(
      X, Y, Cx, Ey, gpParams, numLscAlts, Quants, GPmethod,
      corStruct, optControl, verbose, featuresX, featuresY
) {
    p <- length(featuresX)
    k <- length(featuresY)
    if (missing(gpParams)) {
        if (verbose) {
            message("Fitting GPs for first modality (", p, " features) ...")
        }
        gpsx <- fitManyGPs(
            mat = X, coord = Cx, GPmethod = GPmethod, features = featuresX,
            corStruct = corStruct, optControl = optControl
        )
        if (verbose) {
            message("Fitting GPs for second modality (", k, " features) ...")
        }
        gpsy <- fitManyGPs(
            mat = Y, coord = Ey, GPmethod = GPmethod, features = featuresY,
            corStruct = corStruct, optControl = optControl
        )
    } else {
        # Extract fits
        gpsx <- gpParams$X
        gpsy <- gpParams$Y
    }
    distMat <- as.matrix(stats::dist(rbind(Cx, Ey)))
    n <- nrow(X)
    m <- nrow(Y)
    idN <- seq_len(n)
    idM <- n + seq_len(m) # Indices for x and y
    altSigmas <- buildAltSigmas(distMat,
        numLscAlts = numLscAlts, Quants = Quants,
        idN = idN, idM = idM
    )
    if (verbose) {
        message(
            "Performing all ", p * k,
            " pairwise score tests on fitted GPs ..."
        )
    }
    out <- vapply(selfName(featuresX), function(featx) {
        sx <- base::solve(buildSigmaGp(gpsx[, featx], distMat = distMat[idN, idN]))
        out <- vapply(selfName(featuresY), FUN.VALUE = double(2), function(featy) {
            testGP(
                distMat = distMat, x = X[, featx], y = Y[, featy], altSigmas = altSigmas,
                solXonly = gpsx[, featx], solYonly = gpsy[, featy]
            )
        })
        printProgress(featx, featuresX, verbose)
        return(out)
    }, FUN.VALUE = matrix(0, nrow = 2, ncol = k))
    # Reformat to long format
    t(matrix(c(out), 2, p * k,
        dimnames = list(
            c("pVal", "sign"),
            paste(rep(featuresX, each = k), rep(featuresY, times = p), sep = "__")
        )
    ))
}
