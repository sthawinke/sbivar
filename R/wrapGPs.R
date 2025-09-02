#' Fit Gaussian processes (GPs) if needed, and perform score tests
#'
#' @inheritParams sbivarSingle
#' @inheritParams fitGPs
#' @param gpParams Parameters of the Gaussian processes
#' @returns A named list of results
#' @importFrom smoppix loadBalanceBplapply
#' @importFrom BiocParallel bplapply
wrapGPs = function(X, Y, Cx, Ey, gpParams, families, numLscAlts, Quants, device,
                   GPmethod, training_iter, corStruct, optControl){
    if(missing(gpParams)){
        gpsx = fitManyGPs(mat = X, coord = Cx, family = families[["X"]],
                          device = device, GPmethod = GPmethod,
                          training_iter = training_iter, corStruct = corStruct,
                          optControl = optControl)
        gpsy = fitManyGPs(mat = Y, coord = Ey, family = families[["Y"]],
                          device = device, GPmethod = GPmethod,
                          training_iter = training_iter, corStruct = corStruct,
                          optControl = optControl)
    } else {
        #Extract fits
    }
    distMat = as.matrix(stats::dist(rbind(Cx, Ey)))
    n = nrow(X);m = nrow(Y)
    idN = seq_len(n);idM = n+seq_len(m) #Indices for x and y
    altSigmas = buildAltSigmas(distMat, numLscAlts = numLscAlts, Quants = Quants,
                               idN = idN, idM = idM)
    out = loadBalanceBplapply(selfName(names(gpsx)), function(featx){
        sx = base::solve(buildSigmaGp(gpsx[, featx], distMat = distMat[idN, idN], sparse = FALSE))
        vapply(selfName(names(gpsy)), FUN.VALUE = double(3), function(featy){
            testGP(distMat = distMat, x = X[,featx], y = Y[,featy], altSigmas = altSigmas)
        })
    })
    #Reformat to long format
    t(matrix(unlist(out), 2, length(gpsx)*length(gpsy),
             dimnames = list(c("pVal", "sign"), makeNames(names(gpsx), names(gpsy)))))
}
