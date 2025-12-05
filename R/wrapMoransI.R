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
wrapMoransI = function(X, Y, Cx, Ey, wo, eta, numNN, cutoff, width, verbose, findMaxW, variogramModels, ...){
    if(verbose){
        message("Performing tests on bivariate Moran's I for ", ncol(X)*ncol(Y), " feature pairs")
    }
    #Scale outcomes
    X = scale(X);Y = scale(Y)
    #Move coordinates
    movedCoords = moveTwoCoords(Cx, Ey)
    Cx = movedCoords$Cx;Ey = movedCoords$Ey
    #Test statistic
    W = buildWeightMat(Cx, Ey, wo, eta = eta, numNN = numNN)
    Ixy = (crossprod(X, W) %*% Y)/sqrt(prodFac <- prod(dim(W)-1)) #Normalize for matrix size
    #Variances
    if(verbose){
        message("Fitting variograms for first modality (", ncol(X), " features) ...")
    }
    variogramsX = matheronVariograms(X, Cx, width = width, cutoff = cutoff,
                                     variogramModels = variogramModels,  ...)
    if(verbose){
        message("Fitting variograms for second modality (", ncol(Y), " features) ...")
    }
    variogramsY = matheronVariograms(Y, Ey, width = width, cutoff = cutoff,
                                     variogramModels = variogramModels, ...)
    distX = as.matrix(stats::dist(Cx));distY = as.matrix(stats::dist(Ey))
    if(verbose){
        message("Calculating variances of bivariate Moran's I (", ncol(X)*ncol(Y), " feature pairs) ...")
    }
    varIxy = t(vapply(selfName(colnames(X)), FUN.VALUE = double(ncol(Y)), function(featx){
        sigXw = t(crossprod(W, evalVariogram(variogramsX[[featx]], distX)) %*% W)
        svx = sum(variogramsX[[featx]][, "psill"])
        vapply(selfName(colnames(Y)), FUN.VALUE = double(1), function(featy){
            sum(sigXw * evalVariogram(variogramsY[[featy]], distY))/
                (svx*sum(variogramsY[[featy]][, "psill"]))
            #Fast, memory saving way to find the trace
            #The scaling by variance at the end ensures variances of 1, and is much faster than cov2cor
        })
    }))
    if(any(zeroId <- (varIxy<=0))){
        varIxy[zeroId] = tr(crossprod(W))#If negative variance, fall back on independence
    }
    # P-values
    IxyPvals = makePval(Ixy/sqrt(varIxy/prodFac))
    #Maximum value, if needed
    maxIxy = if(findMaxW) svd(W, nu = 0, nv = 0)$d[1] else 1
    #Reformat to long format
    out <- cbind("Ixy" = c(Ixy), "varIxy" = c(varIxy), "pVal" = c(IxyPvals))
    rownames(out) = makeNames(colnames(X), colnames(Y))
    return(list("out" = out, "maxIxy" = maxIxy))
}
#' Estimate variograms using Matheron's binning estimator for many features at once
#'
#' @importFrom gstat variogram vgm fit.variogram
#' @importFrom sp coordinates
#' @inheritParams sbivarSingle
#' @param X Outcome matrix
#' @param Cx Coordinate matrix
#' @return A list of variograms
#' @details The best fitting variogram model, measured by the squared error, will be used.
matheronVariograms <- function(X, Cx, width, cutoff, variogramModels) {
    df <- data.frame(Cx)
    sp::coordinates(df) <- ~x + y
    # Compute empirical semivariogram using Matheron’s estimator
    variograms <- loadBalanceBplapply(selfName(colnames(X)), function(nm) {
        df$z <- X[, nm]
        vg = variogram(z ~ 1, df, width = width, cutoff = cutoff)
        fvgs = lapply(variogramModels, function(vv){
            fit.variogram(vg, vgm(model = vv, nugget = NA, psill = 1),
                          fit.sills = c(TRUE, FALSE)) #Include nugget variance
        })
        fvg = fvgs[[which.min(vapply(fvgs, FUN.VALUE = double(1), attr, "SSErr"))]]
        if(fvg[2,"range"]<0){
            fvg[2,"range"] = 1e-10 #Catch negative ranges
        }
        return(fvg)
    })
    return(variograms)
}
#' Evaluate a Gaussian variogram on a distance matrix
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
