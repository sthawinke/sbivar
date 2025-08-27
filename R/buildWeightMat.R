#' Build a weight matrix for Moran's I
#'
#' Build different weight matrices to be used in the calculation of the Moran's I
#' statistic.
#'
#' @param xMat,yMat coordinate matrices
#' @param wo The weighting option, either "nn" or "distance
#' @param numNN An integer, the number of neighbours
#' @importFrom RANN nn2
#' @importFrom Matrix sparseMatrix
#' @importFrom spatstat.geom crossdist
#' @return A weight matrix, sparse for the case "nn
buildWeightMat = function(xMat, yMat, wo = c("distance", "nn"), numNN = 8){
    wo = match.arg(wo)
    if(wo == "distance"){
        distMat = spatstat.geom::crossdist(xMat[, 1], xMat[, 2], yMat[, 1], yMat[, 2])
        wm = 1/distMat
        if(ncol(wm)==nrow(wm)){
            diag(wm) = 0
        }
        wm[is.infinite(wm)] = 0
    } else if (wo == "nn"){
        nnMatXY = nn2(xMat, yMat, k = numNN)$nn.idx
        nnMatYX = nn2(yMat, xMat, k = numNN)$nn.idx
        n = nrow(xMat); m = nrow(yMat)
        wm = sparseMatrix(x = rep(1/((n+m)*numNN), (n+m)*numNN),
                          i = c(rep(seq_len(n), times = numNN), c(t(nnMatXY))),
                          j = c(nnMatYX, rep(seq_len(m), each = numNN)))
    }
    return(wm/sum(wm))
}
