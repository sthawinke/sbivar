#' Fit a Gaussian process (GP) to a single outcome vector
#'
#' A Gaussian process (GP) models the spatial autocovariance using a kernel function
#' that describes the covariance as a decreasing function of distance between observations.
#' @note Fitting GPs can be computation and memory intensive!
#' @param x outcome vector
#' @param coord Coordinate matrix
#' @param GPmethod The method by which to fit the Gaussian processes,
#' passed onto \link[nlme]{gls} as "method".
#' @param corStruct The correlation object, see \link[nlme]{corStruct}.
#' At this point, only \link[nlme]{corGaus} is accepted.
#' @param optControl List of control values, see \link[nlme]{glsControl}.
#'
#' @returns A vector of length 4 with components range, nugget, sigma and mean
#' @importFrom stats sigma coef
#' @importFrom nlme gls
fitGP <- function(x, coord, GPmethod, corStruct, optControl) {
    xModGls <- do.call(nlme::gls, args = list(
        "model" = formula("out ~ 1"),
        "data" = data.frame("out" = x, coord),
        "method" = GPmethod, "correlation" = corStruct,
        "control" = optControl
    ))
    xPars <- c(
        coef(xModGls$modelStruct$corStruct, unconstrained = FALSE),
        sigma(xModGls)
    )
    names(xPars) <- c("range", "nugget", "sigma")
    c(xPars, "mean" = unname(coef(xModGls)[1]))
}
#' A wrapper to fit GPs on all columns of a matrix
#'
#' @param mat Matrix of observations
#' @param coord Matrix of coordinates
#' @param ... passed onto \link{fitGP}
#' @inheritParams fitManyGAMs
#' @returns Matrix of fitted GP components
fitManyGPs <- function(mat, coord, features, ...) {
    simplify2array(loadBalanceBplapply(selfName(features), function(cn) {
        fitGP(mat[, cn], coord = coord, ...)
    }))
}
