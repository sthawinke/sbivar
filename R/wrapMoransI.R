#' Calculate ibivariate Moran's I between two modality matrix, with variance and p-value
#'
#' The variance calculation requires estimation of the spatial autocorrelation structure of every feature separately
#'
#' @inheritParams sbivarSingle
#' @param ... passed onto \link[gstat]{variogram}
#'
#' @returns A dataframe of results sorted by p-value, also containing the estimated Moran's I statistic and its variance.
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats dist var
wrapMoransI = function(X, Y, Cx, Ey, wo, eta, numNN, cutoff, width, verbose, ...){
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
    variogramsX = matheronVariograms(X, Cx, width = width, cutoff = cutoff, ...)
    variogramsY = matheronVariograms(Y, Ey, width = width, cutoff = cutoff, ...)
    distX = as.matrix(stats::dist(Cx))
    distY = as.matrix(stats::dist(Ey))
    avDistsX = getAvDists(distX, cutoff, width)
    avDistsY = getAvDists(distY, cutoff, width)
    varIxy = t(simplify2array(loadBalanceBplapply(selfName(colnames(X)), function(featx){
        sigXw = crossprod(W, makeSigmaMatheron(variogramsX[,featx], distId = avDistsX, n = nrow(W))) %*% W
        vapply(selfName(colnames(Y)), FUN.VALUE = double(1), function(featy){
            sigY = makeSigmaMatheron(variogramsY[,featy], distId = avDistsY, n = ncol(W))
            realVar = tr(sigXw %*% sigY)
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
    out <- cbind("Ixy" = c(Ixy), "varIxy" = c(varIxy), "pVal" = c(IxyPvals))
    rownames(out) = makeNames(colnames(X), colnames(Y))
    return(list("out" = out, "maxIxy" = maxIxy))
}
#' Estimate variograms using Matheron's binning estimator for many features at once
#'
#' @importFrom gstat gstat variogram
#' @importFrom sp coordinates
#' @importFrom stats as.formula
#' @param cutoff,width passed onto \link[gstat]{variogram}
matheronVariograms <- function(X, Cx, width, cutoff) {
    df = data.frame(x = Cx[,1], y = Cx[,2], z = X)
    sp::coordinates(df) <- ~x + y
    # Compute empirical semivariogram using Matheron’s estimator
    variograms <- vapply(selfName(colnames(X)), FUN.VALUE = double(floor(cutoff/width)), function(nm) {
        df$z <- X[, nm]
        variogram(z ~ 1, df, width = width, cutoff = cutoff, cross = FALSE)$gamma
    })
    return(variograms)
}
#' Build a spatial autocovariance matrix based on a variogram object
#'
#' @param variogram The variogram, estimated by \link[gstat]{variogram}
#' @param distance The distance matrix
#'
#' @returns The autocovariance matrix
makeSigmaMatheron = function(variogram, distId, n){
    # Convert semivariance to covariance:
    cov_est <- 1 - variogram
    # Interpolate covariance for each distance using nearest distance class
    CovMat <- matrix(cov_est[as.numeric(distId)], nrow = n, ncol = n) #sparseMatrix here?
    diag(CovMat) <- 1
    CovMat[is.na(CovMat)] = 0
    return(CovMat)
}
