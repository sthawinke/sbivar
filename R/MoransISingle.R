#' Calculate bivariate Moran's I between two modality matrix, with variance and p-value
#'
#' The variance calculation requires estimation of the spatial autocorrelation structure of every feature separately
#'
#' @inheritParams sbivarSingle
#' @param ... passed onto \link[gstat]{variogram}
#'
#' @returns A dataframe of results sorted by p-value, also containing the estimated Moran's I statistic and its variance.
#' In addition, the maximum value of the Moran's I statistic, and the parameters of the weight matrix
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats dist
#' @importFrom Matrix forceSymmetric
#'
#' @details By default, a number of range parameters and corresponding weight matrices are screened for spatial association,
#' and their p-value combined using the Cauchy combination rule by \insertCite{Liu2020}{sbivar}.
#'The maximum value of the bivariate Moran's I statistics are returned conditionally,
#' as it is computation intensive and not always needed.
#' @note No multithreading is implemented for the variance calculation, as the matrix calculations involved
#' may use inherent multithreading with OpenBLAS.
MoransISingle = function(X, Y, Cx, Ey, wo, etas, numNN, cutoff, width, verbose,
                         findMaxW, variogramModels, returnSEsMoransI, ...){
    n = nrow(X);m = nrow(Y);p = ncol(X);k=ncol(Y);
    if(verbose)
        message("Testing significance of bivariate Moran's I for ", p*k, " feature pairs")
    #Scale outcomes
    X = scale(X);Y = scale(Y)
    #Move coordinates
    movedCoords = moveTwoCoords(Cx, Ey)
    Cx = movedCoords$Cx;Ey = movedCoords$Ey
    #Estimate spatial autocorrelation
    if(verbose)
        message("Fitting variograms for first modality (", p, " features) ...")
    variogramsX = matheronVariograms(X, Cx, width = width, cutoff = cutoff,
                                     variogramModels = variogramModels,  ...)
    if(verbose)
        message("Fitting variograms for second modality (", k, " features) ...")
    variogramsY = matheronVariograms(Y, Ey, width = width, cutoff = cutoff,
            variogramModels = variogramModels, ...)
    distX = as.vector(stats::dist(Cx));distY = as.vector(stats::dist(Ey))
    prodFac <- (n-1)*(m-1)
    if(verbose)
        message("Calculating bivariate Moran's I statistics ...")
    wParams = selfName(switch(wo, "Gauss" = etas, "nn" = numNN))
    Ws = vapply(wParams, FUN.VALUE = matrix(0, n, m), function(iter) {
        buildWeightMat(Cx = Cx, Ey = Ey, wo = wo, eta = iter, numNN = iter)
    })
    Ws = Ws[,,idW <- (colSums(Ws, dims = 2, na.rm = TRUE)>0), drop = FALSE]
    numWs = dim(Ws)[3]
    if(any(!idW) && (wo=="Gauss")){
        warning("Eta values ", etas[!idW], " yielded zero weight matrices and have been dropped!")
        etas = etas[idW]
    }
    Ixys = vapply(seq_len(numWs), FUN.VALUE = matrix(0, p, k), function(i) {
        crossprod(X, Ws[,,i] %*% Y)
    })/sqrt(prodFac)#Normalize for matrix size
    if(verbose)
        message("Calculating variances of bivariate Moran's I statistics ...")
    #Variances
    mm2 = m*(m-1)/2 #For colSums
    diagMatX = diag(n);ltriX = which(lower.tri(diagMatX))
    ltriY = which(lower.tri(diag(m)))
    varIxy = vapply(selfName(colnames(X)), FUN.VALUE = matrix(0, numWs, k), function(featx){
        diagMatX[ltriX] = evalVariogram(variogramsX[[featx]], distX)
        diagMatX = forceSymmetric(diagMatX, uplo = "L")
        sigXws0 <- lapply(seq_len(numWs), function(i) {
            crossprod(Ws[,,i], diagMatX %*% Ws[,,i])
        }) #BLAS may use multithreading here
        sigXws <- vapply(seq_len(numWs), FUN.VALUE = double(mm2), function(i) {
            sigXws0[[i]][ltriY]
        })
        out = 2*vapply(selfName(colnames(Y)), FUN.VALUE = double(numWs), function(featy){
            vgy = evalVariogram(variogramsY[[featy]], distY)
            .colSums(sigXws*vgy, mm2, numWs) #Fast tr(W^t Sigma_x W Sigma_y)
        }) + vapply(sigXws0, FUN.VALUE = double(1), tr)
        #Diagonal plus two times lower diagonal, exploiting symmetry
        printProgress(featx, colnames(X), verbose)
        return(out)
    })
    varIxy = aperm(varIxy, perm = 3:1) #Rearrange
    for(i in seq_len(numWs)){#If negative variance, fall back on independence
        if(any(zeroId <- (varIxy[,,i]<=0))){
            varIxy[,,i][zeroId] = sum(Ws[,,i]^2) #tr(W^tW)
        }
    }
    varIxy = varIxy/prodFac #Correct for matrix size
    # P-values
    IxyPvals = makePval(Ixys/(seIxy <- sqrt(varIxy)))
    #CCT correction
    cctPvals = apply(IxyPvals, c(1,2), CCT)
    #Reformat to long format
    out <- matrix(c(Ixys), ncol = numWs, dimnames = list(NULL, paste0("Ixy_", wParams)))
    if(returnSEsMoransI){
      out = cbind(out, matrix(c(seIxy), ncol = numWs,
                    dimnames = list(NULL, paste0("SE(Ixy)_", wParams))))
    }
    out = cbind(out, "pVal" = c(cctPvals))
    rownames(out) = makeNames(colnames(X), colnames(Y))
    #Maximum values, if needed
    maxIxy = if(findMaxW) {
        vapply(selfName(names(wParams)), FUN.VALUE = double(1), function(i) {
            svd(Ws[,,i], nu = 0, nv = 0)$d[1]
        })
    }
    return(list("res" = out, "wo" = wo, "etas" = if(wo=="Gauss") wParams,
                "maxIxy" = maxIxy, "numNN" = if(wo=="nn") wParams))
}
#' Estimate variograms using Matheron's binning estimator for many features at once, and evaluate
#'
#' @importFrom gstat variogram vgm fit.variogram
#' @importFrom sp coordinates
#' @inheritParams sbivarSingle
#' @param X Outcome matrix
#' @param Cx Coordinate matrix
#' @return A list of evaluated variograms
#' @details The best fitting variogram model, measured by the squared error, will be used.
matheronVariograms <- function(X, Cx, width, cutoff, variogramModels) {
    Cx = data.frame(Cx)
    sp::coordinates(Cx) <- ~x + y
    # Compute empirical semivariogram using Matheron’s estimator
    variograms <- loadBalanceBplapply(selfName(colnames(X)), function(nm) {
        Cx$z <- X[, nm]
        vg = variogram(z ~ 1, Cx, width = width, cutoff = cutoff)
        fvgs = lapply(variogramModels, function(vv){
            try(fit.variogram(vg, vgm(model = vv, nugget = NA)), silent = TRUE) #Include nugget variance
        })
        fvgs = fvgs[vapply(fvgs, FUN.VALUE = TRUE, is, "variogramModel")]
        fvg = fvgs[[which.min(vapply(fvgs, FUN.VALUE = double(1), attr, "SSErr"))]]
        if(fvg[2,"range"]<0){
            fvg[2,"range"] = 1e-10 #Catch negative ranges
            fvg[2,"psill"] = 0
        }
        fvg[,"psill"] = fvg[,"psill"]/sum(fvg[,"psill"]) #Normalize to variance 1
        return(fvg)
    })
    return(variograms)
}
#' Evaluate a variogram on a set of covariances
#'
#' @param vg The variogram
#' @param distVec A vector of distances
#' @returns A vector of covariances
evalVariogram = function(vg, distVec){
    covMat = vg[2, "psill"]*if(vg[2, "model"] == "Exp"){
        exp(-distVec/vg[2, "range"])
    } else if(vg[2, "model"] == "Lin"){
        tmp = numeric(length(distVec))
        id = distVec < vg[2, "range"]
        tmp[id] = 1-distVec[id]/vg[2, "range"]
        tmp
    }
    return(covMat)
}


