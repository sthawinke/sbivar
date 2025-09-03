#' Fit a gaussian process to a single outcome vector
#'
#' @param x outcome vector
#' @param coord Coordinate matrix
#' @param device Device to fit the GPs on. Defaults to "cpu",
#' can use "gpu" via the python package gpytorch if installed
#' @param GPmethod The method by which to fit the Gaussian processes,
#' passed onto \link[nlme]{gls} as "method" if it equals "REML" or "ML".
#' @param training_iter Number of training iterations in gpytorch
#' @param corStruct The correlation object, see \link[nlme]{corStruct}.
#' At this point, only \link[nlme]{corGaus} is accepted
#' @param optControl List of control values, see \link[nlme]{glsControl}
#' @details Providing "cuda" as device exploits the available GPU for speedup,
#' especially for large GPs (thousands of observations). Setting GPmethod to
#'  "ML" or "REML" leads to model fitting in R using the \link[nlme]{gls} function.
#' Setting GPmethod to "gpytorch" uses the eponymous python package, with the option of gpu acceleration.
#' The parametrization of the GPs is the one common in geostatistics, as described in \link[nlme]{gls}.
#'
#' @returns A vector of length 4 with components range, nugget, sigma and mean
#' @importFrom stats sigma coef
#' @importFrom nlme gls
fitGPs = function(x, coord, GPmethod, device, training_iter, corStruct, optControl){
    if(GPmethod == "gpytorch"){
        pyFit = fitGPsPython(x, coord, training_iter = training_iter, device = device)
        return(unlist(pyFit$Pars))
    } else {
        xModGls <- do.call(nlme::gls, args = list("model" = formula("out ~ 1"),
                            "data" = data.frame("out" = x, coord),
                            "method" = GPmethod, "correlation" = corStruct,
                            "control" = optControl))
        xPars = c(coef(xModGls$modelStruct$corStruct, unconstrained = FALSE),
                  sigma(xModGls))
        names(xPars) = c("range", "nugget", "sigma")
        c(xPars, "mean" = unname(coef(xModGls)[1]))
    }
}
#' A wrapper to fit GPs on all columns of a matrix
#'
#' @param mat Matrix of observations
#' @param coord Matrix of coordinates
#' @param ... passed onto \link{fitGPs}
#' @returns Matrix of fitted GP components
fitManyGPs = function(mat, coord, ...){
    simplify2array(loadBalanceBplapply(selfName(colnames(mat)), function(cn){
        fitGPs(mat[, cn], coord = coord, ...)
    }))
}
