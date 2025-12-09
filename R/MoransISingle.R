#' Calculate ibivariate Moran's I between two modality matrix, with variance and p-value
#'
#' The variance calculation requires estimation of the spatial autocorrelation structure of every feature separately
#'
#' @inheritParams sbivarSingle
#' @param ... passed onto \link[gstat]{variogram}
#'
#' @returns A dataframe of results sorted by p-value, also containing the estimated Moran's I statistic and its variance.
#' In addition, the maximum value of the Moran's I statistic
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats dist
#'
#' @details The maximum value of the bivariate Moran's I statistic is returned conditionally,
#' as it is computation intensive and not always needed.
MoransISingle = function(X, Y, Cx, Ey, wo, etas, numNN, cutoff, width, verbose, findMaxW, variogramModels, ...){
    n = nrow(X);m = nrow(Y);p = ncol(X);k=ncol(Y);
    if(verbose){
        message("Testing significance of bivariate Moran's I for ", p*k, " feature pairs")
    }
    #Scale outcomes
    X = scale(X);Y = scale(Y)
    #Move coordinates
    movedCoords = moveTwoCoords(Cx, Ey)
    Cx = movedCoords$Cx;Ey = movedCoords$Ey
    #Estimate spatial autocorrelation
    if(verbose){
        message("Fitting variograms for first modality (", p, " features) ...")
    }
    variogramsX = matheronVariograms(X, Cx, width = width, cutoff = cutoff,
                                     variogramModels = variogramModels,  ...)
    if(verbose){
        message("Fitting variograms for second modality (", k, " features) ...")
    }
    variogramsY = matheronVariograms(Y, Ey, width = width, cutoff = cutoff,
            variogramModels = variogramModels, ...)
    distX = as.matrix(stats::dist(Cx));distY = as.matrix(stats::dist(Ey))
    prodFac <- (n-1)*(m-1)
    #Weight matrices and test statistics
    Ws = vapply(switch(wo, "Gauss" = etas, "nn" = numNN), FUN.VALUE = matrix(0, n, m), function(iter) {
        buildWeightMat(Cx = Cx, Ey = Ey, wo = wo, eta = iter, numNN = iter)
    })
    Ws = Ws[,,idW <- (colSums(Ws, dims = 2)>0), drop = FALSE]
    numWs = dim(Ws)[3]
    if(any(!idW) && wo=="Gauss"){
        warning("Eta values ", etas[!idW], " yielded zero weight matrices and have been dropped!")
        etas = etas[idW]
    }
    Ixys = vapply(seq_len(numWs), FUN.VALUE = matrix(0, p, k), function(i) {
        (crossprod(X, Ws[,,i]) %*% Y)/sqrt(prodFac) #Normalize for matrix size
    })
    if(verbose){
        message("Calculating variances of bivariate Moran's I (", p*k, " feature pairs) ...")
    }
    #Variances
    ncs = m^2 #For colSums
    varIxy = vapply(selfName(colnames(X)), FUN.VALUE = matrix(0, numWs, k), function(featx){
        vgx = evalVariogram(variogramsX[[featx]], distX)
        sigXws = vapply(seq_len(numWs), FUN.VALUE = matrix(0, m, m), function(i) {
            crossprod(Ws[,,i], vgx) %*% Ws[,,i]
        }) #BLAS may use multithreading here
        out = simplify2array(loadBalanceBplapply(selfName(colnames(Y)), function(featy){
            .colSums(sigXws*c(evalVariogram(variogramsY[[featy]], distY)), ncs, numWs)
            #Fast, memory saving way to find the trace
        }))
        if(verbose)
            printProgress(featx, colnames(X))
        return(out)
    })
    varIxy = aperm(varIxy, perm = 3:1) #Rearrange
    for(i in seq_len(numWs)){
        if(any(zeroId <- (varIxy[,,i]<=0))){
            varIxy[,,i][zeroId] = tr(crossprod(Ws[,,i]))
            #If negative variance, fall back on independence
        }
    }
    # P-values
    IxyPvals = makePval(Ixys/sqrt(varIxy/prodFac))
    #CCT correction
    cctPvals = apply(IxyPvals, c(1,2), CCT)
    #Reformat to long format
    out <- cbind(matrix(c(Ixys), ncol = numWs, dimnames = list(NULL, paste0("Ixy_",
                switch(wo, "Gauss" = etas, "nn" = numNN)))), "pVal" = c(cctPvals))
    rownames(out) = makeNames(colnames(X), colnames(Y))
    #Maximum values, if needed
    maxIxy = if(findMaxW) {
        vapply(seq_len(numWs), FUN.VALUE = double(1), function(i) {
            svd(Ws[,,i], nu = 0, nv = 0)$d[1]
        })
    }
    return(list("out" = out, "etas" = if(wo=="Gauss") etas, "maxIxy" = maxIxy,
                "numNN" = if(wo=="nn") numNN))
}
#' Estimate variograms using Matheron's binning estimator for many features at once, and evaluate
#'
#' @importFrom gstat variogram vgm fit.variogram
#' @importFrom sp coordinates
#' @inheritParams sbivarSingle
#' @param X Outcome matrix
#' @param Cx Coordinate matrix
#' @return An array of evaluated variograms
#' @details The best fitting variogram model, measured by the squared error, will be used.
matheronVariograms <- function(X, Cx, width, cutoff, variogramModels) {
    Cx = data.frame(Cx)
    sp::coordinates(Cx) <- ~x + y
    # Compute empirical semivariogram using Matheron’s estimator
    variograms <- loadBalanceBplapply(selfName(colnames(X)), function(nm) {
        Cx$z <- X[, nm]
        vg = variogram(z ~ 1, Cx, width = width, cutoff = cutoff)
        fvgs = lapply(variogramModels, function(vv){
            fit.variogram(vg, vgm(model = vv, nugget = NA)) #Include nugget variance
        })
        fvg = fvgs[[which.min(vapply(fvgs, FUN.VALUE = double(1), attr, "SSErr"))]]
        if(fvg[2,"range"]<0){
            fvg[2,"range"] = 1e-10 #Catch negative ranges
            fvg[2,"psill"] = 0
        }
        fvg[,"psill"] = fvg[,"psill"]/sum(fvg[,"psill"]) #Normalize to variance 1
        return(data.frame(fvg[, c("model", "psill", "range")]))
    })
    return(variograms)
}
#' Evaluate a variogram on a distance matrix
#'
#' @param vg The variogram
#' @param distMat The distance matrix
#' @returns A covariance matrix
evalVariogram = function(vg, distMat){
    covMat = vg[2, "psill"]*if(vg[2, "model"] == "Gau"){
        exp(-(distMat/vg[2, "range"])^2)
    } else if(vg[2, "model"] == "Sph"){
        dvg = distMat/vg[2, "range"]
        tmp =(1-1.5*dvg+0.5*dvg^3)
        tmp[distMat > vg[2, "range"]] = 0
        tmp
    } else if(vg[2, "model"] == "Exp"){
        exp(-distMat/vg[2, "range"])
    } else if(vg[2, "model"] == "Lin"){
        tmp = 1-distMat/vg[2, "range"]
        tmp[distMat > vg[2, "range"]] = 0
        tmp
    }
    diag(covMat) = diag(covMat) + vg[1, "psill"]
    return(covMat)
}
