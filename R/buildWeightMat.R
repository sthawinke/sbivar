#' Build a weight matrix for bivariate Moran's I
#'
#' Build a weight matrix to be used in the calculation of the bivariate Moran's I
#' statistic.
#'
#' @param Cx,Ey Two coordinate matrices
#' @param wo The weighting option, either "nn" or "Gauss"
#' @param numNN An integer, the number of neighbours
#' @param eta parameter that controls the decay of the weights with distance, see details
#' @param distMat The distance matrix
#' @importFrom Matrix sparseMatrix
#' @importFrom spatstat.geom crossdist
#' @return A weight matrix
#' @details
#' For wo = "Gauss", the weight decays as exp(-d^2/eta) with d the distance between observations.
#' For wo = "nn" the numNN nearest neighbours are given equal weight, all others are zero.
buildWeightMat <- function(
      Cx, Ey, wo, eta, numNN,
      distMat = spatstat.geom::crossdist(Cx[, 1], Cx[, 2], Ey[, 1], Ey[, 2])
) {
    if (wo == "nn") {
        if (requireNamespace("RANN")) {
            nnMatXY <- RANN::nn2(Cx, Ey, k = numNN)$nn.idx
            nnMatYX <- RANN::nn2(Ey, Cx, k = numNN)$nn.idx
            n <- nrow(Cx)
            m <- nrow(Ey)
            wm <- as.matrix(sparseMatrix(
                x = rep(1 / ((n + m) * numNN), (n + m) * numNN),
                i = c(rep(seq_len(n), times = numNN), c(t(nnMatXY))),
                j = c(nnMatYX, rep(seq_len(m), each = numNN))
            ))
        } else {
            stop("Install 'RANN' package for nearest neighbour weighting matrices")
        }
    } else if (wo == "Gauss") {
        wm <- exp(-distMat^2 / eta)
    }
    return(wm / sum(wm))
}
