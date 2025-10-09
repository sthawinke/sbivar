#' Find all Moran's I statistics for a list of matrices
#' @description
#' Find all raw Moran's I values between lists of observations matrices from different modalities.
#' Observation matrices are scaled and centered prior to calculation
#'
#' @inheritParams sbivarMulti
#'
#' @returns A list of named Moran's I vectors
#' @seealso \link{buildWeightMat}
wrapMoransIMulti = function(Xl, Yl, Cxl, Eyl, wo, numNN, verbose){
    lapply(selfName(names(Xl)), function(nam){
        if(verbose)
            printIteration(nam, names(Xl))
        out = c(crossprod(scale(Xl[[nam]]), buildWeightMat(Cxl[[nam]], Eyl[[nam]], wo = wo, numNN = numNN))
                %*% scale(Yl[[nam]]))
        names(out) = makeNames(colnames(Xl[[nam]]), colnames(Yl[[nam]]))
        return(out)
    })
}
