#' Title
#'
#' @inheritParams sbivarSingle
#' @param x outcome vector
#' @param GPmethod The method by which to fit the Gaussian processes
#' @param training_iter Number of training iterations in gpytorch
#' @param device Device on which the GPs are fitted using the python package \textit{gpytorch}
#' @param corStruct The correlation object, see \link[nlme]{corStruct}
#' @param optControl List of control values, see \link[nlme]{glsControl}
#' @details Providing "cuda" as device exploits the available GPU for speedup,
#' especially for large GPs (thousands of observations). Setting GPmethod to
#'  "ML" or "REML" leads to model fitting in R using the \link[nlme]{gls} function.
#' Setting GPmethod to "gpytorch" uses the eponymous python package, with the option of gpu accceleration.
#' The parametrization of the GPs is the one common in geostatistics, as described in \link[nlme]{gls}.
#'
#' @returns A vector of length 4 with components range, nugget, sigma and mean
#' @importFrom stats sigma coef
#' @importFrom nlme gls
fitGPs = function(x, Cx, GPmethod, training_iter, device, corStruct, optControl){
    if(GPmethod == "gpytorch"){
        pyFit = fitGPsPython(xMat, x, training_iter = training_iter, device = device)
        return(unlist(pyFit$Pars))
    } else {
        xModGls <- do.call(nlme::gls, args = list("model" = formula("out ~ 1"), "data" = data.frame("out" = x, Cx),
                                                  "method" = GPmethod, "correlation" = corStruct, "control" = optControl))
        xPars = c(coef(xModGls$modelStruct$corStruct, unconstrained = FALSE), sigma(xModGls))
        names(xPars) = c("range", "nugget", "sigma")
        c(xPars, "mean" = coef(xModGls)[1])
    }
}
#' A wrapper to fit GPs on all columns of a matrix
#'
#' @inheritParams sbivarSingle
#' @param ... passed onto \link{fitGPs}
#'
#' @returns Matrix of fitted GP components
fitManyGPs = function(X, Cx, ...){
    simplify2array(loadBalanceBplapply(selfName(colnames(X)), function(cn){
        fitGPs(X[, cn], Cx = Cx, ...)
    }))
}
