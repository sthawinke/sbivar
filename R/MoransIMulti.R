#' Find all Moran's I statistics for a list of matrices
#' @description
#' Find all raw Moran's I values between lists of observations matrices from different modalities.
#' Observation matrices are scaled and centered prior to calculation
#'
#' @inheritParams sbivarMulti
#' @inheritParams buildWeightMat
#' @param ... passed onto \link{MoransISingle}
#'
#' @returns A list of named Moran's I vectors
#' @seealso \link{buildWeightMat}
MoransIMulti = function(Xl, Yl, Cxl, Eyl, verbose, ...){
    lapply(selfName(names(Xl)), function(nam){
        if(verbose)
            printIteration(nam, names(Xl))
        MoransISingle(X = Xl[[nam]], Y = Yl[[nam]], Cx = Cxl[[nam]], Ey = Eyl[[nam]],
              verbose = FALSE, findMaxW = TRUE, returnVarsMoransI = TRUE, ...)
    })
}
