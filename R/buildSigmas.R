#' Construct cross-blocks of alternative covariance matrices at different length scales
#'
#' These matrices are used to test for bivariate association at different length scales.
#' Each alternative covariance matrix has the block structure
#' \eqn{\Sigma_{alt,l} = \begin{bmatrix} I_n & C_l \\ C_l^T & I_m \end{bmatrix}}
#' where \eqn{C_l} is the Gaussian cross-covariance kernel evaluated at length scale \eqn{l}.
#' Only the \eqn{n \times m} cross-blocks are returned; the identity diagonal is implicit.
#'
#' @inheritParams testGP
#' @param numLscAlts Number of length scales (and thus number of cross-blocks to be built)
#' @param Quants Most extreme quantiles of the distance distribution to be used as length scales.
#' @param idN,idM indices for x and y in the distance matrix
#'
#' @returns An \eqn{n \times m \times L} array of cross-blocks \eqn{C_l}
#' @importFrom stats quantile
buildAltSigmas <- function(distMat, numLscAlts, Quants, idN, idM) {
    rangeDist <- {
        tmp <- distMat[upper.tri(distMat)]
        quantile(tmp[tmp != 0], probs = Quants)
    }
    lscAlts <- exp(seq(log(rangeDist[1]), log(rangeDist[2]), length.out = numLscAlts))
    crossDist <- distMat[idN, idM]   # n x m sub-matrix of distances
    # Return only the off-diagonal cross-blocks C_l = GaussKernel(crossDist, lscAlt_l)
    # Result is an n x m x L array; the diagonal identity blocks are implicit.
    vapply(lscAlts, FUN.VALUE = crossDist, function(lscAlt) {
        GaussKernel(crossDist, lscAlt)
    })
}
#' Build the SAC matrix for a Gaussian process
#'
#' @param pars A vector of parameters
#' @inheritParams testGP
#'
#' @returns The covariance matrix
#' @seealso \link[nlme]{corGaus}, \link[nlme]{corMatrix}
buildSigmaGp <- function(pars, distMat) {
    pars["sigma"]^2 * (diag(nrow(distMat)) * pars["nugget"] +
        (1 - pars["nugget"]) * GaussKernel(distMat, pars["range"]))
}
#' Construct the SAC part of the covariance matrix using the Gaussian covariance kernel,
#'
#' @inheritParams testGP
#' @param range Range (or length-scale) parameter of the Gaussian covariance kernel
#'
#' @returns The SAC covariance matrix
GaussKernel <- function(distMat, range) {
    exp(-(distMat / range)^2)
}
