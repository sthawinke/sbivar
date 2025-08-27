#' Return predictions of a gam, along with their variance-covariance matrix
#'
#' The covariance matrix of the predictions is found using the formula:
#' Cov(pred) = X * Var(coef) * X'
#' Where X is the basis matrix (design matrix for smooth terms), and Var(coef) is the covariance matrix of the coefficients
#' From help(predict.gam): When type="lpmatrix" then a matrix is returned which
#' yields the values of the linear predictor (minus any offset) when postmultiplied by
#' the parameter vector (in this case se.fit is ignored). The latter option is most useful
#' for getting variance estimates for quantities derived from the model: for example integrated quantities, or derivatives of smooths.
#'
#' @param model The gam model
#' @param newdata The grid on which new predictions are done
#' @return A list with components
#' \item{pred}{A vector of predictions}
#' \item{vcov}{Variance-covariance matrix of the predictions}
#' @seealso \link[mgcv]{vcov.gam}, \link[mgcv]{predict.gam}
#' @importFrom mgcv vcov.gam predict.gam
vcovPredGam = function(model, newdata){
    # Get the full coefficient covariance matrix
    coef_cov_matrix <- vcov.gam(model)
    # Get the basis matrix for predictions
    basis_matrix <- predict.gam(model, newdata = newdata, type = "lpmatrix")
    predOut = c(basis_matrix %*% coef(model))
    #Same as predict(model, newdata = newdata, type = "response")
    predOut = switch(model$family$link,
                     "identity" = predOut, "inverse" = 1/predOut, "log" = exp(predOut))
    prediction_cov_matrix <- basis_matrix %*% tcrossprod(coef_cov_matrix, basis_matrix)
    #Checked using predict.gam(,se.fit = TRUE)
    return(list("pred" = predOut, "vcov" = prediction_cov_matrix))
}
#' Get approximate variance-covariance matrix of spline predictions
#'
#' @param vcovMat Variance-covariance matrix of spline coefficients
#' @param cen Centered observations
#' @param x Observations
#' @param link Link function
#'
#' @returns Variance-covariance matrix of the spline predictions
getApproxVar = function(vcovMat, cen, x, link){
    vec = switch(link, "identity" = cen, "log" = cen*x, "inverse" = cen/x^2)
    return(vec %*% vcovMat %*% vec)
    #Minus for inverse link (gamma) cancels out
}
