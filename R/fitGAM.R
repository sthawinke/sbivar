#' Fit a GAM model to a single variable
#'
#' A generalized additive model (GAM) captures spatial dependence through a bivariate spatial spline.
#' @param df The dataframe containing outcome and coordinates
#' @param outcome A character vector indicating the outcome variable
#' @param k Dimension of the basis of the smooth term see \link[mgcv]{s}.
#' @param family A character string indicating the family, see \link[stats]{family}.
#' @param offset A numeric vector with the offset, e.g. the library sizes for
#' spatial transcriptomic data, already on the scale of the regressor,
#' so log-transformed for count models
#' @returns A fitted GAM model, or try-error when the fit fails
#' @importFrom mgcv gam s
#' @import stats
#' @seealso \link[mgcv]{gam}
fitGAM = function(df, outcome, k = -1, family = gaussian(), offset = NULL){
    try(gam(as.formula(paste(outcome, " ~ s(x, y, k = k)")), data = df, family = family,
        offset = offset), silent = TRUE)
}
#' Fit GAMs to all columns of a dataframe, as a wrapper for fitGAM
#'
#' @param mat The matrix of outcomes
#' @param coord The coordinate matrix
#' @param modality Character vector indicating which modality is being fit.
#' For debugging purposes mainly
#' @param ... Passed onto \link{fitGAM}
#' @inheritParams fitGAM
#'
#' @returns A list of GAM models
#' @importFrom smoppix loadBalanceBplapply
#' @importFrom BiocParallel bplapply
fitManyGAMs = function(mat, coord, family = gaussian(), modality, ...){
    cns = selfName(colnames(mat))
    df = data.frame(as.matrix(mat), coord)
    if(family$family != "gaussian"){
        libSizes = rowSums(mat)
        df = df[id <- (libSizes > 0),]
        libSizes = libSizes[id]
    }
    offset = switch(family$link, "inverse" = 1/libSizes,
                    "log" = log(libSizes), NULL)
    fits <- loadBalanceBplapply(cns, function(cn){
        fitGAM(df, outcome = cn, offset = offset, family = family, ...)
    })
    if(!any(id <- vapply(fits, FUN.VALUE = TRUE, is, "gam"))){
        stop("All GAM fits failed in modality ", modality,
             ", please investigate cause! First error message:\n", fits[[1]])
    }
    return(fits[id])
}

