#' Perform nearest-neighbour matching, and apply all pairwise tests
#'
#' @inheritParams sbivar
#' @param mapToFinest A boolean, should the mapping occur to the dataset with
#' the best resolution?
#'
#' @returns A dataframe of results sorted by p-value
#' @importFrom RANN nn2
#' @importFrom SpatialPack modified.ttest
wrapModTtest = function(X, Y, Cx, Ey, mapToFinest = FALSE){
    n = nrow(X);m = nrow(Y)
    if(xor(n > m, mapToFinest)){
        idnn_XY <- nn2(Cx, Ey, k = 1)$nn.idx[, 1]
        X = X[idnn_XY,]
        coordMat = Cx[idnn_XY,]
    } else {
        idnn_YX <- nn2(Ey, Cx, k = 1)$nn.idx[, 1]
        Y = Y[idnn_YX,];
        coordMat = Ey[idnn_YX,]
    }
    featGrid = expand.grid("featX" = colnames(X), "featY" = colnames(Y))
    out = simplify2array(loadBalanceBplapply(seq_len(nrow(featGrid)), function(i){
        modified.ttest(X[, featGrid[i, "featX"]], Y[, featGrid[i, "featY"]],
                       coordMat)[c("corr", "dof", "p.value")]
    }))
    colnames(out) = apply(featGrid, 1, paste, collapse = "_")
    data.frame(t(out))
}
