#' Fit Gaussian processes (GPs) if needed, and perform score tests
#'
#' @inheritParams sbivarSingle
#' @param gpParams Parameters of the Gaussian processes
wrapGPs = function(X, Y, Cx, Ey, gpParams, families, numLscAlts, Quants, ...){
    xGPs = fitManyGPs(mat = X, coord = Cx, family = families[["X"]], ...)
    yGPs = fitManyGPs(mat = Y, coord = Ey, family = families[["Y"]], ...)
    testManyGPs(xGPs, yGPs, numLscAlts = numLscAlts, Quants = Quants)
}
