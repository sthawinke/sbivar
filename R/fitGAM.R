#' Fit a GAM model to a single variable
#'
#' Fit a generalized additive model (GAM) captures spatial dependence through a bivariate spatial spline.
#' @param df The dataframe containing outcome and coordinates
#' @param outcome A character vector indicating the outcome variable
#' @param k Dimension of the basis of the smooth term, see \link[mgcv]{s}.
#' @param family A character string indicating the family, see \link[stats]{family}.
#' @param offset A numeric vector with the offset, e.g. the library sizes for
#' spatial transcriptomic data, already on the scale of the regressor,
#' so log-transformed for count models.
#' @returns A fitted GAM model, or try-error when the fit fails
#' @importFrom mgcv gam gamm s
#' @import stats
#' @details If a gamma fit is attempted and fails, which frequenlty happens for sparse data,
#' a negative binomial fit is attempted instead
#' @seealso \link[mgcv]{gam},\link[mgcv]{s}
#' @inheritParams fitGP
fitGAM <- function(df, outcome, k = -1, family = gaussian(), offset = NULL, Gamm, correlation) {
    Form <- as.formula(paste(outcome, " ~ s(x, y, k = k)"))
    fit <- if(Gamm){
        try(gamm(Form, correlation = correlation,
                data = df, family = family,
                offset = offset)$gam, silent = TRUE)
    } else {
        try(gam(Form,
        data = df, family = family,
        offset = offset), silent = TRUE)
    }
    if (is(fit, "try-error") && family$family == "Gamma") {
        fit <- fitGAM(
            df = df, outcome = outcome, k = k, family = mgcv::nb(),
            offset = offset
        )
    }
    return(fit)
}
#' Fit GAMs to all columns of a dataframe, as a wrapper for fitGAM
#'
#' @param mat The matrix of outcomes
#' @param coord The coordinate matrix
#' @param modality Character vector indicating which modality is being fit.
#' For debugging purposes mainly
#' @param ... Passed onto \link{fitGAM}
#' @inheritParams fitGAM
#' @param pseudoCount Pseudocount added to avoid zeroes for gamma distribution
#' @param features Features to be fit, others are only used to estimate the offset
#'
#' @returns A list of GAM models
#' @importFrom smoppix loadBalanceBplapply
#' @importFrom BiocParallel bplapply
fitManyGAMs <- function(mat, coord, family = gaussian(), modality, features, pseudoCount = 1e-8, ...) {
    if (family$family == "Gamma") {
        mat <- mat + pseudoCount
    }
    df <- data.frame(as.matrix(mat), coord)
    if (family$family != "gaussian") {
        libSizes <- rowSums(mat)
        df <- df[id <- (libSizes > 0), ]
        libSizes <- libSizes[id]
    }
    offset <- switch(family$link,
        "inverse" = 1 / libSizes,
        "log" = log(libSizes),
        NULL
    )
    fits <- loadBalanceBplapply(selfName(features), function(cn) {
        fitGAM(df, outcome = cn, offset = offset, family = family, ...)
    })
    if (!all(id <- vapply(fits, FUN.VALUE = TRUE, is, "gam"))) {
        warning(
            immediate. = TRUE,
            sum(!id), " GAM fits failed in modality ", modality,
            ", please investigate cause! First failure:\n", fits[[which.min(id)]]
        )
    }
    return(fits[id])
}
