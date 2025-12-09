#' Find all cross-correlations for a list of matrices
#' @description
#' Find all raw cross-correlations between lists of observations matrices from different modalities.
#'
#' @inheritParams sbivarMulti
#' @inheritParams wrapModTtest
#'
#' @returns A list of named correlation vectors
#' @importFrom stats cor
correlationsMulti = function(Xl, Yl, verbose){
    lapply(selfName(names(Xl)), function(nam){
        if(verbose)
            printIteration(nam, names(Xl))
        out = c(cor(Xl[[nam]], Yl[[nam]]))
        names(out) = makeNames(colnames(Xl[[nam]]), colnames(Yl[[nam]]))
        return(out)
    })
}
