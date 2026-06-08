#' Find all cross-correlations for a list of matrices
#' @description
#' Find all raw cross-correlations between lists of observations matrices from different modalities.
#' @param verbose Should progress be printed?
#' @inheritParams sbivarMulti
#'
#' @returns A list of named correlation vectors
#' @importFrom stats cor
correlationsMulti <- function(Xl, Yl, featuresX, featuresY, verbose) {
    lapply(selfName(names(Xl)), function(nam) {
        if (verbose) {
            printIteration(nam, names(Xl))
        }
        commonNames <- intersect(rownames(Xl[[nam]]), rownames(Yl[[nam]]))
        # Ensure same dimensions
        out <- c(cor(
            Xl[[nam]][commonNames, featuresX <- intersect(featuresX, colnames(Xl[[nam]]))],
            Yl[[nam]][commonNames, featuresY <- intersect(featuresY, colnames(Yl[[nam]]))]
        ))
        names(out) <- makeNames(featuresX, featuresY)
        return(list("res" = out))
    })
}
