#' Test for correlation between the predictions of two GAM models
#'
#' @param modelx,modely Two fitted GAMs
#' @param predx,predy Predictions and covariance matrices of fitted GAMs in common grid
#'
#' @returns A vector with correlation, its standard error and the p-value
testGAM = function(modelx, modely, predx, predy){
    cxy = sum((cen1 <- (predx$pred-mean(predx$pred)))*(cen2 <- (predy$pred-mean(predy$pred))))
    #Covariance
    approxVar = getApproxVar(predx$vcov, cen2, predx$pred, modelx$family$link) +
        getApproxVar(predy$vcov, cen1, predy$pred, modely$family$link)
    #Its variance, exploit block diagonality
    se = sqrt(approxVar)
    denom = (length(predx$pred)-1)*sd(predx$pred)*sd(predy$pred)
    corxy = cxy/denom
    #Correlation
    return(c("corxy" = corxy, "se.corxy" = se/denom, "pVal" = makePval(cxy/se)))
}
#' Test all bivariate combinations for fitted lists of GAMs
#'
#' @param gamsx,gamsy Lists of GAM models for two modalities
#' @param newGrid The new grid of points in which to evaluate the GAMs
#'
#' @returns A named list of results
#' @importFrom smoppix loadBalanceBplapply
#' @importFrom BiocParallel bplapply
testManyGAMs = function(gamsx, gamsy, newGrid){
    gamsx = gamsx[!vapply(gamsx, FUN.VALUE = TRUE, inherits, "try-error")]
    gamsy = gamsy[!vapply(gamsy, FUN.VALUE = TRUE, inherits, "try-error")]
    loadBalanceBplapply(selfName(names(gamsx)), function(featx){
        predx <- vcovPredGam(gamsx[[featx]], newdata = newGrid)
        vapply(selfName(names(gamsy)), FUN.VALUE = double(3), function(featy){
            testGAM(predx = predx, modely = gamsy[[featy]], modelx = gamsx[[featx]],
                         predy = vcovPredGam(gamsy[[featy]], newdata = newGrid))
        })
    })
}
#' Wrapper function to fit GAMs and test for all possible combinations
#'
#' @inheritParams sbivar
#'
#' @returns see \link{testManyGAMs}
wrapGAMs = function(X, Y, Cx, Ey, families, n_points_grid){
    xGams = fitManyGAMs(mat = X, coord = Cx, family = families[["X"]])
    yGams = fitManyGAMs(mat = Y, coord = Ey, family = families[["Y"]])
    ng = buildNewGrid(Cx = Cx, Ey = Ey, n_points_grid = n_points_grid)
    testManyGAMs(xGams, yGams, newGrid = ng)
}
