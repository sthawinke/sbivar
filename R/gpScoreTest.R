#' Perform the score test for significance of the bivariate spatial association for Gaussian processes.
#'
#' @inheritParams sbivarSingle
#' @param x,y outcome vectors
#' @param device Device to fit the GPs on. Defaults to "cpu",
#' can use "gpu" via the python package gpytorch if installed
#' @param altSigmas A prepared series of bivariate association matrices
#' @param method passed onto \link[nlme]{gls}
#' @param distMat The distance matrix of Cx and Ey
#' @param solXonly,solYonly Parametersof the Gaussian process
#' @param sx Inverse of the covariance matrix of x. Will be calculated if missing.
#' @inheritParams buildAltSigmas
#'
#' @returns A vector of length 2: a p-value and an indicator of the sign:
#' +1 for positive association, -1 for negative
#' @importFrom stats pchisq
#' @importFrom abind abind
gpScoreTest = function(x, y, Cx, Ey,
                       device = "cpu", altSigmas, method = "gaussian", distMat = as.matrix(stats::dist(rbind(Cx, Ey))),
                       solXonly, solYonly, sx, numLscAlts = if(missing(altSigmas)) 10 else dim(altSigmas)[3],
                       Quants = c(0.005, 0.5)){
    n = length(x);m=length(y)
    #Estimate separate GPs
    if(missing(solXonly))
        solXonly = prepGLS(x, Cx, method = method, device = device)
    if(missing(solYonly))
        solYonly = prepGLS(y, Ey, method = method, device = device)
    out = c(x, y); muVec = rep(c(solXonly["mean"], solYonly["mean"]), times = c(n,m))
    diffVec = cbind(out-muVec)
    idN = seq_len(n);idM = n+seq_len(m) #Indices for x and y
    #Exploit block diagonality
    if(missing(sx))
        sx <- base::solve(buildSigmaGp(solXonly, distMat =  distMat[idN, idN], sparse = FALSE))
    invW <- bdiag(sx, base::solve(buildSigmaGp(solYonly, distMat = distMat[idM, idM], sparse = FALSE)))
    P = invW - (invW %*% regMat %*% solve(crossprod(regMat, invW) %*% regMat, crossprod(regMat,  invW)))
    #Need derivatives for the efficient information
    derivList = list("X" = abind("sigmaX" = arrayDeriv(solXonly, distMat[idN, idN], what = "sigmax"),
                                 "rangeX" = arrayDeriv(solXonly, distMat[idN, idN], what = "range"),
                                 "nuggetX" = arrayDeriv(solXonly, distMat[idN, idN], what = "nugget"), along = 3),
                     "Y" = abind("sigmaY" = arrayDeriv(solYonly, distMat[idM, idM], what = "sigmax"),
                                 "rangeY" = arrayDeriv(solYonly, distMat[idM, idM], what = "range"),
                                 "nuggetY" = arrayDeriv(solYonly, distMat[idM, idM], what = "nugget"), along = 3))
    #Make two arrays out of this list
    derivListP = list("X" = arrayMatProd(M = P[idN, idN], A = derivList[["X"]]),
                      "Y" = arrayMatProd(M = P[idM, idM], A = derivList[["Y"]]))
    Ithetatheta = 0.5*bdiagn(arrayProd(derivListP$X, derivListP$X), arrayProd(derivListP$Y, derivListP$Y))
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
    if(missing(altSigmas))
        altSigmas = buildAltSigmas(distMat, numLscAlts = numLscAlts, Quants = Quants, idN = idN, idM = idM)
    e = 0.5*tr(P)
    pSigma = arrayMatProd(M = P, A = altSigmas)
    Itautau = 0.5*colSums(pSigma^2, dims = 2)
    pSigma = arrayMatProd(M = P, A = pSigma)
    Itautheta = 0.5*t(rbind(arrayProd2tr(A = pSigma[idN, idN,], B = derivList$X),
                            arrayProd2tr(A = pSigma[idM, idM,], B = derivList$Y)))
    ItautauTilde = Itautau - rowSums((Itautheta[,idItt] %*% sitt)* Itautheta[,idItt])
    kappaEst = ItautauTilde/(2*e);nu = 2*e^2/ItautauTilde
    pVals = vapply(c(FALSE, TRUE), FUN.VALUE = double(numLscAlts), function(neg){
        #Regular variances goes into the weights (invW), only the alternative sigma is in the test statistic
        ## The score statistic
        vec = c(if(neg) invWDiffVecNeg else invWDiffVec)
        Usigma = 0.5 * crossprod(colSums(vec * altSigmas), vec)
        pchisq(as.vector(Usigma/kappaEst), df = nu, lower.tail = FALSE)
    })
    #Cauchy combination rule
    ps = apply(pVals, 2, function(p) CCT(p[!is.na(p)]))
    idMin = which.min(ps)
    c("ScorePval" = min(2*ps[idMin], 1), "sign" = if(idMin==1) 1 else -1)
    #two one-sided tests, multiply p-value by 2 to get two-sided test
}
