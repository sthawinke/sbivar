#' Fit Gaussian processes (GPs) if needed, and perform score tests
#'
#' @inheritParams sbivarSingle
#' @param gpParams Parameters of the Gaussian processes
wrapGPs = function(X, Y, Cx, Ey, gpParams, families, numLscAlts, Quants, device,
                   GPmethod, training_iter, corStruct, optControl){
    if(missing(gpParams)){
        xGPs = fitManyGPs(mat = X, coord = Cx, family = families[["X"]],
                          device = device, GPmethod = GPmethod,
                          training_iter = training_iter, corStruct = corStruct,
                          optControl = optControl)
        yGPs = fitManyGPs(mat = Y, coord = Ey, family = families[["Y"]],
                          device = device, GPmethod = GPmethod,
                          training_iter = training_iter, corStruct = corStruct,
                          optControl = optControl)
    }
    testManyGPs(xGPs, yGPs, numLscAlts = numLscAlts, Quants = Quants, X = X, Y = Y,
                distMat = as.matrix(stats::dist(rbind(Cx, Ey))))
}
