#' Find all cross-correlations for a list of matrices
#' @description
#' Find all raw cross-correlations between lists of observations matrices from different modalities.
#' Observations are matched to nearest neighbours if necessary
#'
#' @inheritParams sbivarMulti
#' @inheritParams wrapModTtest
#'
#' @returns A list of named correlation vectors
#' @importFrom stats cor
#' @seealso \link{matchCoords}
wrapCorrelationsMulti = function(Xl, Yl, Cxl, Eyl, mapToFinest, jointCoordinates){
    lapply(selfName(names(Xl)), function(nam){
        if(!jointCoordinates && !identical(Xl[[nam]], Yl[[nam]])){
            matchedCoords = matchCoords(Cxl[[nam]], Eyl[[nam]], mapToFinest = mapToFinest)
            coordMat = matchedCoords$coordMat
            if(matchedCoords$mapToX){
                Xl[[nam]] = Xl[[nam]][matchedCoords$id,]
            } else {
                Yl[[nam]] = Yl[[nam]][matchedCoords$id,]
            }
        } else {
            coordMat = Cxl[[nam]]
        }
        out = c(cor(Xl[[nam]], Yl[[nam]]))
        names(out) = makeNames(colnames(Xl[[nam]]), colnames(Yl[[nam]]))
        return(out)
    })
}
