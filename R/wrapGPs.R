#' Fit Gaussian processes (GPs) if needed, and perform score tests
#'
#' @inheritParams sbivarSingle
#' @inheritParams fitGP
#' @param numLscAlts Number of length scales to be tested for bivariate association
#' @param Quants Most extreme quantiles of the distance distribution to take as length scales
#' @returns A named list of results
#' @importFrom smoppix loadBalanceBplapply
#' @importFrom BiocParallel bplapply
wrapGPs = function(X, Y, Cx, Ey, gpParams, numLscAlts, Quants, device,
                   GPmethod, training_iter, corStruct, optControl){
    if(missing(gpParams)){
        if(GPmethod == "gpytorch"){
            stop("Interface with python not yet implemented, set GP method to 'REML' or 'ML'")
        }
        gpsx = fitManyGPs(mat = X, coord = Cx,
                          device = device, GPmethod = GPmethod,
                          training_iter = training_iter, corStruct = corStruct,
                          optControl = optControl)
        gpsy = fitManyGPs(mat = Y, coord = Ey,
                          device = device, GPmethod = GPmethod,
                          training_iter = training_iter, corStruct = corStruct,
                          optControl = optControl)
    } else {
        #Extract fits
        gpsx = gpParams$X;gpsy = gpParams$Y
    }
    distMat = as.matrix(stats::dist(rbind(Cx, Ey)))
    n = nrow(X);m = nrow(Y)
    idN = seq_len(n);idM = n+seq_len(m) #Indices for x and y
    altSigmas = buildAltSigmas(distMat, numLscAlts = numLscAlts, Quants = Quants,
                               idN = idN, idM = idM)
    out = loadBalanceBplapply(selfName(colnames(X)), function(featx){
        sx = base::solve(buildSigmaGp(gpsx[, featx], distMat = distMat[idN, idN]))
        vapply(selfName(colnames(Y)), FUN.VALUE = double(2), function(featy){
            testGP(distMat = distMat, x = X[,featx], y = Y[,featy], altSigmas = altSigmas,
                   solXonly = gpsx[, featx], solYonly = gpsy[, featy])
        })
    })
    #Reformat to long format
    t(matrix(unlist(out), 2, ncol(X)*ncol(Y),
             dimnames = list(c("pVal", "sign"),
    paste(rep(colnames(X), each = ncol(Y)), rep(colnames(Y), times = ncol(X)), sep = "__"))))
}
