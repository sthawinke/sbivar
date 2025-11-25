#' Calculate ibivariate Moran's I between two modality matrix, with variance and p-value
#'
#' The variance calculation requires estimation of the spatial autocorrelation structure of every feature separately
#'
#' @inheritParams sbivarSingle
#'
#' @returns A dataframe of results sorted by p-value, also containing the estimated Moran's I statistic and its variance.
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats dist var
wrapMoransI = function(X, Y, Cx, Ey, wo, eta, numNN, verbose, ...){
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
    Ixy = (crossprod(X, W) %*% Y)/sqrt(prod(dim(W)-1)) #Normalize for matrix size
    #Variances
    variogramsX = matheronVariograms(X, Cx, ...)
    variogramsY = matheronVariograms(Y, Ey, ...)
    distX = as.matrix(stats::dist(Cx))
    distY = as.matrix(stats::dist(Ey))
    varsIxy = loadBalanceBplapply(selfName(colnames(X)), function(featx){
        sigXw = crossprod(W, makeSigmaMatheron(variogramsX[variogramsX$id == featx,], distance = distX)) %*% W
        vapply(selfName(colnames(Y)), FUN.VALUE = double(1), function(featy){
            sigY = makeSigmaMatheron(variogramsY[variogramsY$id == featy,], distance = distY)
            realVar = tr(sigXw %*% sigY)
            return(realVar)
        })
    })
    indepVar = tr(crossprod(W))#If negative variance, fall back on independence
    varsIxy[varsIxy<=0] = indepVar
    # P-values
    IxyPvals = makePval(Ixy/sqrt(varsIxy))
    #Maximum value
    maxIxy = svd(W, nu = 0, nv = 0)$d[1]
    #Reformat to long format
    out = cbind("Ixy" = c(Ixy), "varIxy" = c(varIxy), "pVal" = c(IxyPvals))
    rownames(out) = makeNames(colnames(X), colnames(Y))
    attr(out, "max") = maxIxy
    return(out)
}
#' @importFrom gstat gstat variogram
#' @importFrom sp coordinates
#' @param cutoff,width passed onto \link[gstat]{variogram}
matheronVariograms <- function(X, Cx, cutoff, width) {
    df = data.frame(x = coords[,1], y = coords[,2], z = X)
    sp::coordinates(df) <- ~x + y
    g <- NULL
    for (nm in colnames(X)) {
        g <- gstat(g, id = nm, formula = as.formula(paste0(nm,"~1")), data=df)
    }
    # Compute empirical semivariogram using Matheron’s estimator
    v <- variogram(g, cutoff = cutoff, width = width, cross = FALSE)[,c("dist", "gamma", "id")]
    return(v)
}
#' Build a spatial autocovariance matrix based on a variogram object
#'
#' @param variogram The variogram, estimated by \link[gstat]{variogram}
#' @param distance The distance matrix
#'
#' @returns The autocovariance matrix
makeSigmaMatheron = function(variogram, distance){
    # Convert semivariance to covariance:
    cov_est <- 1 - variogram$gamma
    # Interpolate covariance for each distance using nearest distance class
    nearest_index <- findInterval(distance, variogram$dist)
    # Avoid 0 index
    nearest_index[nearest_index == 0] <- 1
    CovMat <- matrix(cov_est[nearest_index], nrow = nrow(coords), ncol = nrow(coords))
    diag(CovMat) <- 1
    return(CovMat)
}
