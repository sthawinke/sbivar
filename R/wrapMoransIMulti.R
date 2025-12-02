#' Find all Moran's I statistics for a list of matrices
#' @description
#' Find all raw Moran's I values between lists of observations matrices from different modalities.
#' Observation matrices are scaled and centered prior to calculation
#'
#' @inheritParams sbivarMulti
#' @inheritParams buildWeightMat
#'
#' @returns A list of named Moran's I vectors
#' @seealso \link{buildWeightMat}
wrapMoransIMulti = function(Xl, Yl, Cxl, Eyl, wo, numNN, eta, verbose){
    lapply(selfName(names(Xl)), function(nam){
        if(verbose)
            printIteration(nam, names(Xl))
        #Bring coordinates to 0-1
        minX = min(c(Cxl[[nam]][, "x"], Eyl[, "x"]))
        minY = min(c(Cxl[[nam]][, "y"], Eyl[, "y"]))
        Cxl[[nam]][, "x"] = Cxl[[nam]][, "x"] - minX
        Cxl[[nam]][, "y"] = Cxl[[nam]][, "y"] - minY
        Eyl[[nam]][, "x"] = Eyl[[nam]][, "x"] - minX
        Eyl[[nam]][, "y"] = Eyl[[nam]][, "y"] - minY
        MaxCoord = max(c(Cxl[[nam]], Eyl[[nam]]))
        Cxl[[nam]] = Cxl[[nam]]/MaxCoord;Eyl[[nam]] = Eyl[[nam]]/MaxCoord
        out = c(crossprod(scale(Xl[[nam]]), buildWeightMat(Cxl[[nam]], Eyl[[nam]],
                eta = eta, wo = wo, numNN = numNN)) %*% scale(Yl[[nam]]))
        names(out) = makeNames(colnames(Xl[[nam]]), colnames(Yl[[nam]]))
        return(out)
    })
}
