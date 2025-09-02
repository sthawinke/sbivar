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

