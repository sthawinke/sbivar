#' Perform nearest-neighbour matching if necessary, and apply all pairwise tests.
#'
#' If measurements are not on the same location, they are matched using nearest neighbour matching with the \link[RANN]{nn2} function.
#' Then modified t-test is applied to all pairs
#'
#' @inheritParams sbivarSingle
#' @param jointCoordinates A boolean, are measurements on the same location
#'.
#' @returns A dataframe of results sorted by p-value, also containing effective sample size (ESS) and correlation estimate.
#' @importFrom RANN nn2
#' @importFrom SpatialPack modified.ttest
#' @seealso \link[SpatialPack]{modified.ttest}
wrapModTtest = function(X, Y, Cx, Ey, mapToFinest = FALSE, jointCoordinates = FALSE){
    n = nrow(X);m = nrow(Y)
    if(!jointCoordinates){
        if(xor(n > m, mapToFinest)){
            idnn_XY <- nn2(Cx, Ey, k = 1)$nn.idx[, 1]
            X = X[idnn_XY,]
            coordMat = Cx[idnn_XY,]
        } else {
            idnn_YX <- nn2(Ey, Cx, k = 1)$nn.idx[, 1]
            Y = Y[idnn_YX,];
            coordMat = Ey[idnn_YX,]
        }
    } else {
        coordMat = Cx
    }
    featGrid = expand.grid("featX" = colnames(X), "featY" = colnames(Y))
    out = simplify2array(loadBalanceBplapply(seq_len(nrow(featGrid)), function(i){
        unlist(modified.ttest(X[, featGrid[i, "featX"]], Y[, featGrid[i, "featY"]],
                       coordMat)[c("corr", "ESS", "p.value")])
    }))
    colnames(out) = apply(featGrid, 1, paste, collapse = "_")
    rownames(out) = c("Correlation", "ESS", "pVal")
    t(out)
}
