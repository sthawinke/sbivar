#' Perform the score test for significance of the bivariate spatial association for Gaussian processes.
#'
#' @inheritParams sbivarSingle
#' @param x,y outcome vectors
#' @param device Device to fit the GPs on. Defaults to "cpu",
#' can use "gpu" via the python package gpytorch if installed
#' @param altSigmas A prepared series of bivariate association matrices
#' @param GPmethod passed onto \link[nlme]{gls} as "method"
#' @param distMat The distance matrix of Cx and Ey
#' @param solXonly,solYonly Parametersof the Gaussian process
#' @param sx Inverse of the covariance matrix of x. Will be calculated if missing.
#' @inheritParams buildAltSigmas
#'
#' @returns A vector of length 2: a p-value and an indicator of the sign:
#' +1 for positive association, -1 for negative
#' @importFrom stats pchisq dist
#' @importFrom abind abind
#' @importFrom Matrix bdiag
testGP = function(x, y, Cx, Ey, device = "cpu", altSigmas, GPmethod,
                       distMat = as.matrix(stats::dist(rbind(Cx, Ey))),
                       solXonly, solYonly, sx, numLscAlts = if(missing(altSigmas)) 10 else dim(altSigmas)[3],
                       Quants = c(0.005, 0.5)){
    n = length(x);m=length(y)
    idN = seq_len(n);idM = n+seq_len(m) #Indices for x and y
    regMat = cbind(1, c(rep(0, n), rep(1, m))) #Regression matrix accounting for mean differences
    #Estimate separate GPs
    if(missing(solXonly))
        solXonly = fitGPs(x, Cx, GPmethod = GPmethod, device = device)
    if(missing(solYonly))
        solYonly = fitGPs(y, Ey, GPmethod = GPmethod, device = device)
    out = c(x, y); muVec = rep(c(solXonly["mean"], solYonly["mean"]), times = c(n,m))
    diffVec = cbind(out-muVec)

    #Exploit block diagonality
    if(missing(sx))
        sx <- base::solve(buildSigmaGp(solXonly, distMat =  distMat[idN, idN], sparse = FALSE))
    invW <- bdiag(sx, base::solve(buildSigmaGp(solYonly, distMat = distMat[idM, idM], sparse = FALSE)))
    P = invW - (invW %*% regMat %*% solve(crossprod(regMat, invW) %*% regMat, crossprod(regMat,  invW)))
    #Need derivatives for the efficient information
    derivList = list("X" = abind("sigmaX" = arrayDeriv(solXonly, distMat[idN, idN], what = "sigma"),
                                 "rangeX" = arrayDeriv(solXonly, distMat[idN, idN], what = "range"),
                                 "nuggetX" = arrayDeriv(solXonly, distMat[idN, idN], what = "nugget"), along = 3),
                     "Y" = abind("sigmaY" = arrayDeriv(solYonly, distMat[idM, idM], what = "sigma"),
                                 "rangeY" = arrayDeriv(solYonly, distMat[idM, idM], what = "range"),
                                 "nuggetY" = arrayDeriv(solYonly, distMat[idM, idM], what = "nugget"), along = 3))
    #Make two arrays out of this list
    derivListP = list("X" = arrayMatProd(M = P[idN, idN], A = derivList[["X"]]),
                      "Y" = arrayMatProd(M = P[idM, idM], A = derivList[["Y"]]))
    Ithetatheta = 0.5*bdiagn(arrayProd(derivListP$X, derivListP$X),
                             arrayProd(derivListP$Y, derivListP$Y))
    sitt = try(solve(Ithetatheta), silent = TRUE)
    idItt = colnames(sitt)
    if(inherits(sitt, "try-error")){
        idItt = grep("range", colnames(Ithetatheta), invert = TRUE)
        sitt = try(solve(Ithetatheta[idItt, idItt]), silent = TRUE)
    } #If singular, ignore range parameter,
    if(inherits(sitt, "try-error")){
        idItt = grep("sigma", colnames(Ithetatheta))
        sitt = try(solve(Ithetatheta[idItt, idItt]), silent = TRUE)
    } #If singular, ignore range and sigma parameters
    invWDiffVec = invWDiffVecNeg = invW %*% diffVec #Prepare matrix
    invWDiffVecNeg[idM] = -invWDiffVecNeg[idM] #Prepare matrix
    e = 0.5*tr(P)
    pSigma = arrayMatProd(M = P, A = altSigmas)
    Itautau = 0.5*colSums(pSigma^2, dims = 2)
    pSigma = arrayMatProd(M = P, A = pSigma)
    Itautheta = 0.5*t(rbind(arrayProd2tr(A = pSigma[idN, idN,], B = derivList$X),
                            arrayProd2tr(A = pSigma[idM, idM,], B = derivList$Y)))
    ItautauTilde = Itautau - rowSums((Itautheta[,idItt] %*% sitt)* Itautheta[,idItt])
    kappaEst = ItautauTilde/(2*e);nu = 2*e^2/ItautauTilde
    pVals = vapply(c(FALSE, TRUE), FUN.VALUE = double(numLscAlts), function(neg){
        vec = c(if(neg) invWDiffVecNeg else invWDiffVec)
        Usigma = 0.5 * crossprod(colSums(vec * altSigmas), vec) ## The score statistic
        pchisq(as.vector(Usigma/kappaEst), df = nu, lower.tail = FALSE)
    })
    #Cauchy combination rule
    ps = apply(pVals, 2, function(p) CCT(p[!is.na(p)]))
    idMin = which.min(ps)
    c("pVal" = min(2*ps[idMin], 1), "sign" = if(idMin==1) 1 else -1)
    #two one-sided tests, multiply p-value by 2 to get two-sided test
}
#' Test all bivariate combinations for fitted lists of GAMs
#'
#' @param gpsx,gpssy Lists of GAM models for two modalities
#' @param newGrid The new grid of points in which to evaluate the GAMs
#'
#' @returns A named list of results
#' @importFrom smoppix loadBalanceBplapply
#' @importFrom BiocParallel bplapply
testManyGPs = function(gpsx, gpsy, X, Y, distMat){
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
#' Construct the nxnxg/2 array of derivatives for a nxn matrix to the g/2 covariance matrix parameters.
#'
#' @param fittedGP The fitted Gaussian process (vector of 4 parameters)
#' @inheritParams gpScoreTest
#' @param what For which parameter is the derivative required?
#' @return The matrix of derivatives
arrayDeriv = function(fittedGP, distMat, what){
    switch(what,
           "sigma" = 2*fittedGP["sigma"]*(exp(-(distMat/fittedGP["range"])^2)*(1-fittedGP["nugget"]) +
                                              diag(fittedGP["nugget"], nrow(distMat))),
           "nugget" = {
               tmp = -exp(-(distMat/fittedGP["range"])^2)*fittedGP["sigma"]^2
               diag(tmp) = 0
               tmp
           },
           "range" = {
               tmp = 2*(1-fittedGP["nugget"])*fittedGP["sigma"]^2*
                   exp(-(distMat/fittedGP["range"])^2)*distMat^2/fittedGP["range"]^3
               diag(tmp) = 0
               tmp
           })
}
