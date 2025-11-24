#' Build a weight matrix for bivariate Moran's I
#'
#' Build a weight matrix to be used in the calculation of the bivariate Moran's I
#' statistic.
#'
#' @param Cx,Ey Two coordinate matrices
#' @param wo The weighting option, either "nn" or "distance"
#' @param numNN An integer, the number of neighbours
#' @param eta parameter that controls the decay of the weights with distance, see details
#' @importFrom RANN nn2
#' @importFrom Matrix sparseMatrix
#' @importFrom spatstat.geom crossdist
#' @return A weight matrix
#' @details
#' For wo = "exp", the weight decays as exp(-d^2/eta) with d the distance between observations.
buildWeightMat = function(Cx, Ey, wo = c("exp", "distance", "nn"), eta = 0.05, numNN = 8){
    wo = match.arg(wo)
    if(wo == "distance"){
        distMat = spatstat.geom::crossdist(Cx[, 1], Cx[, 2], Ey[, 1], Ey[, 2])
        wm = 1/distMat
        if(ncol(wm)==nrow(wm)){
            diag(wm) = 0
        }
        wm[is.infinite(wm)] = 0
    } else if (wo == "nn"){
        nnMatXY = nn2(Cx, Ey, k = numNN)$nn.idx
        nnMatYX = nn2(Ey, Cx, k = numNN)$nn.idx
        n = nrow(Cx); m = nrow(Ey)
        wm = as.matrix(sparseMatrix(x = rep(1/((n+m)*numNN), (n+m)*numNN),
                          i = c(rep(seq_len(n), times = numNN), c(t(nnMatXY))),
                          j = c(nnMatYX, rep(seq_len(m), each = numNN))))
    } else if (wo == "exp"){
        distMat = spatstat.geom::crossdist(xMat[, 1], xMat[, 2], yMat[, 1], yMat[, 2])
        wm = exp(-distMat^2/eta)
    }
    return(wm/sum(wm))
}
