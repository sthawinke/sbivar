#' Title
#'
#' @inheritParams gpScoreTest
#' @param numLscAlts Number of length scales (and thus number of covariance matrices to be tested)
#' @param Quants Most extreme quantiles of the distance distribution to be used as length scales.
#' @param idN,idM indices for x and y in the distance matrix
#'
#' @returns A list of covariance matrices
#' @importFrom stats quantile
#' @importFrom Matrix sparseMatrix
buildAltSigmas = function(distMat, numLscAlts, Quants, idN, idM){
    rangeDist = {
        tmp = distMat[upper.tri(distMat)]
        quantile(tmp[tmp!=0], probs = Quants)
    }
    lscAlts = exp(seq(log(rangeDist[1]), log(rangeDist[2]), length.out = numLscAlts))
    n = length(idN);m = length(idM)
    vapply(lscAlts, FUN.VALUE = distMat, function(lscAlt) {
        #No sigma needed here, just the skeleton
        mat = sparseMatrix(i = rep(idN, times = m), j = rep(idM, each = n), dims = c(n+m, n+m),
                           symmetric = TRUE, x = c(exp(-(distMat[idN, idM]/lscAlt)^2)))
        diag(mat) <- 1
        return(as.matrix(mat))
    })
}

