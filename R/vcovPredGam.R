#' Return predictions of a GAM, along with the factored coefficient covariance
#'
#' Spline predictions are \eqn{B \beta}, where \eqn{B} is the basis (design) matrix
#' and \eqn{\beta} are the spline coefficients.  Rather than materialising the full
#' \eqn{N_\text{grid} \times N_\text{grid}} prediction covariance
#' \eqn{B \, \text{Var}(\beta) \, B^T}, only \eqn{B} (\eqn{N_\text{grid} \times q})
#' and \eqn{\text{Var}(\beta)} (\eqn{q \times q}) are returned.  The quadratic form
#' \eqn{v^T B \, \text{Var}(\beta) B^T v = (B^T v)^T \text{Var}(\beta) (B^T v)}
#' is then computed cheaply in \eqn{q} dimensions by \link{getApproxVar}.
#'
#' @param model  The fitted GAM
#' @param newdata The grid on which predictions are made
#' @return A list with components
#' \item{pred}{A vector of predictions on the response scale}
#' \item{basis}{\eqn{N_\text{grid} \times q} basis matrix for the smooth of interest}
#' \item{coef_cov}{\eqn{q \times q} unconditional covariance matrix of the
#'   smooth coefficients}
#' @seealso \link[mgcv]{vcov.gam}, \link[mgcv]{predict.gam}
#' @importFrom mgcv vcov.gam predict.gam
#' @inheritParams MoransISingle
#' @inheritParams GAMsSingle
vcovPredGam <- function(model, newdata, findVariances = TRUE) {
    # Identify the smooth of interest and its coefficient indices
    idSmooth <- which(vapply(model$smooth,
        FUN.VALUE = TRUE,
        function(sm) sm$id == "trend"
    ))
    idCoefs <- seq(
        model$smooth[[idSmooth]]$first.para,
        model$smooth[[idSmooth]]$last.para
    )

    # Basis matrix B (N_grid x q) — shared by prediction and variance computation
    basis_matrix <- predict.gam(model,
        newdata = newdata, type = "lpmatrix",
        newdata.guaranteed = TRUE
    )[, idCoefs]

    predOut <- c(basis_matrix %*% coef(model)[idCoefs])
    predOut <- switch(model$family$link,
        "identity" = predOut,
        "inverse"  = 1 / predOut,
        "log"      = exp(predOut)
    )

    if (findVariances) {
        # q x q covariance of smooth coefficients — much smaller than N_grid x N_grid
        coef_cov <- vcov.gam(model, unconditional = TRUE)[idCoefs, idCoefs]
        return(list("pred" = predOut, "basis" = basis_matrix, "coef_cov" = coef_cov))
    }
    return(list("pred" = predOut))
}

#' Approximate variance of the spatial covariance numerator
#'
#' Computes \eqn{v^T \text{Cov}(\hat{f}) \, v} using the factored form
#' \eqn{(B^T v)^T C (B^T v)}, where \eqn{B} is the basis matrix and
#' \eqn{C = \text{Var}(\beta)}. This avoids materialising the
#' \eqn{N_\text{grid} \times N_\text{grid}} prediction covariance matrix.
#'
#' @param predInfo List returned by \link{vcovPredGam} (must contain
#'   \code{basis} and \code{coef_cov})
#' @param cen  Centered predictions of the \emph{other} modality
#' @param x    Raw predictions (needed for non-identity links)
#' @param link Link function name
#'
#' @returns Scalar approximate variance
getApproxVar <- function(predInfo, cen, x, link) {
    vec <- switch(link,
        "identity" = cen,
        "log"      = cen * x,
        "inverse"  = cen / x^2
    )
    # u = B^T vec  (q x 1, q << N_grid)
    u <- crossprod(predInfo$basis, vec)
    # u^T C u  (scalar quadratic form in q dimensions)
    c(crossprod(u, predInfo$coef_cov %*% u))
}
