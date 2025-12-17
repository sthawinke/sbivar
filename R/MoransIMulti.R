#' Find all Moran's I statistics for a list of matrices
#'
#' Find all Moran's I values for lists of observations matrices from different modalities,
#' by caling the \link{MoransISingle} function.
#'
#' @inheritParams sbivarMulti
#' @inheritParams buildWeightMat
#' @param ... passed onto \link{MoransISingle}
#'
#' @returns A list of Moran's I estimates, standard errors and maximum values
#' @seealso \link{MoransISingle}
MoransIMulti = function(Xl, Yl, Cxl, Eyl, returnSEsMoransI, verbose, ...){
    lapply(selfName(names(Xl)), function(nam){
        if(verbose)
            printIteration(nam, names(Xl))
        MoransISingle(X = Xl[[nam]], Y = Yl[[nam]], Cx = Cxl[[nam]], Ey = Eyl[[nam]],
              verbose = FALSE, findMaxW = TRUE, returnSEsMoransI = returnSEsMoransI, ...)[c("res", "maxIxy")]
    })
}
