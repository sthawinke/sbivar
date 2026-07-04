#' Fit a GAM model to a single variable
#'
#' Fit a generalized additive model (GAM) captures spatial dependence through a bivariate spatial spline.
#' @param df The dataframe containing outcome and coordinates
#' @param outcome A character vector indicating the outcome variable
#' @param family A character string indicating the family, see \link[stats]{family}.
#' @param offset A numeric vector with the offset, e.g. the library sizes for
#' spatial transcriptomic data, already on the scale of the regressor,
#' so log-transformed for count models.
#' @returns A fitted GAM model, or try-error when the fit fails
#' @importFrom mgcv gam gamm s
#' @import stats
#' @details If a gamma fit is attempted and fails, which frequently happens for sparse data,
#' a negative binomial fit is attempted instead
#' @seealso \link[mgcv]{gam}, \link[mgcv]{s}
fitGAM <- function(df, outcome, family = gaussian(), Gamm, correlation, offset = NULL) {
    Form <- paste(outcome, "~ s(x, y, bs = 'tp', id = 'trend')",
                  if(Gamm && !is.null(df$Offset)) "+offset(Offset)")
    fit <- if(Gamm){
        try(gamm(as.formula(Form), correlation = correlation,
                data = df, family = family
        )$gam, silent = TRUE)
        } else {
            try(gam(as.formula(Form),
        data = df, family = family, offset = df$Offset,
    ), silent = TRUE)
        }
    if (is(fit, "try-error") && family$family == "Gamma") {
        fit <- fitGAM(
            df = df, outcome = outcome, family = mgcv::nb(),
            Gamm = Gamm, correlation = correlation
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
fitManyGAMs <- function(mat, coord, family = gaussian(), modality, features,
                        Gamm, correlation, pseudoCount = 1e-8, ...) {
    if (family$family == "Gamma") {
        mat <- mat + pseudoCount
    }
    df <- data.frame(as.matrix(mat), coord)
    if (family$family != "gaussian") {
        libSizes <- rowSums(mat)
        df <- df[id <- (libSizes > 0), ]
        libSizes <- libSizes[id]
    }
    df$Offset <- switch(family$link,
        "inverse" = 1 / libSizes,
        "log" = log(libSizes),
        NULL
    )
    fits <- loadBalanceBplapply(selfName(features), function(cn) {
        fitGAM(df, outcome = cn, family = family, Gamm = Gamm, correlation = correlation, ...)
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
