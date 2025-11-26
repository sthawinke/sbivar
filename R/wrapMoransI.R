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
#' @importFrom gstat variogramLine
#'
#' @details The Moran's I values and variances returned have been scaled by their maximum value,
#' the first singular vector of the weight matrix., already.
wrapMoransI = function(X, Y, Cx, Ey, wo, eta, numNN, cutoff, width, verbose, model, ...){
    if(verbose){
        message("Performing tests on bivariate Moran's I for ", ncol(X)*ncol(Y), " feature pairs")
    }
    #Scale outcomes
    X = scale(X);Y = scale(Y)
    #Bring coordinates to 0-1
    minX = min(c(Cx[, "x"], Ey[, "x"]));minY = min(c(Cx[, "y"], Ey[, "y"]))
    Cx[, "x"] = Cx[, "x"] - minX; Cx[, "y"] = Cx[, "y"] - minY
    Ey[, "x"] = Ey[, "x"] - minX; Ey[, "y"] = Ey[, "y"] - minY
    MaxCoord = max(c(Cx, Ey))
    Cx = Cx/MaxCoord;Ey = Ey/MaxCoord
    #Test statistic
    W = buildWeightMat(Cx, Ey, wo, eta = eta, numNN = numNN)
    Ixy = (crossprod(X, W) %*% Y)/sqrt(prodFac <- prod(dim(W)-1)) #Normalize for matrix size
    #Variances
    variogramsX = matheronVariograms(X, Cx, width = width, cutoff = cutoff, model = model, ...)
    variogramsY = matheronVariograms(Y, Ey, width = width, cutoff = cutoff, model = model, ...)
    distX = as.matrix(stats::dist(Cx))
    distY = as.matrix(stats::dist(Ey))
    varIxy = t(simplify2array(loadBalanceBplapply(selfName(colnames(X)), function(featx){
        sigXw = crossprod(W, variogramLine(variogramsX[[featx]], dist_vector = distX, covariance = TRUE)) %*% W
        vapply(selfName(colnames(Y)), FUN.VALUE = double(1), function(featy){
            realVar = tr(sigXw %*% variogramLine(variogramsY[[featy]], dist_vector = distY, covariance = TRUE))
            return(realVar)
        })
    })))
    indepVar = tr(crossprod(W))#If negative variance, fall back on independence
    varIxy[varIxy<=0] = indepVar
    # P-values
    IxyPvals = makePval(Ixy/sqrt(varIxy/prodFac))
    #Maximum value
    maxIxy = svd(W, nu = 0, nv = 0)$d[1]
    #Reformat to long format
    out <- cbind("Ixy" = c(Ixy/maxIxy), "varIxy" = c(varIxy/maxIxy), "pVal" = c(IxyPvals))
    rownames(out) = makeNames(colnames(X), colnames(Y))
    return(list("out" = out, "maxIxy" = maxIxy))
}
#' Estimate variograms using Matheron's binning estimator for many features at once
#'
#' @importFrom gstat variogram vgm fit.variogram
#' @importFrom sp coordinates
#' @inheritParams sbivarSingle
matheronVariograms <- function(X, Cx, width, cutoff, model) {
    df = data.frame(x = Cx[,1], y = Cx[,2], z = X)
    sp::coordinates(df) <- ~x + y
    # Compute empirical semivariogram using Matheron’s estimator
    variograms <- loadBalanceBplapply(selfName(colnames(X)), function(nm) {
        df$z <- X[, nm]
        fvg = fit.variogram(variogram(z ~ 1, df, width = width, cutoff = cutoff),
                      vgm(1, model = model), fit.sills = FALSE) #Fix partial sill at 1, includes nugget variance
        if(fvg$range<0){
            fvg$range = 1e-10 #Catch negative ranges
        }
        return(fvg)
    })
    return(variograms)
}
