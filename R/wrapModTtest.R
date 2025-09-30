#' Perform nearest-neighbour matching if necessary, and apply all pairwise tests.
#'
#' If measurements are not on the same location, they are matched using nearest neighbour matching with the \link[RANN]{nn2} function.
#' Then modified t-test is applied to all pairs, which tests for the significance
#' of the Pearson correlation while accounting for spatial autocorrelation \insertCite{Clifford1989}{sbivar}.
#'
#' @inheritParams sbivarSingle
#' @param jointCoordinates A boolean, are measurements on the same location
#'.
#' @returns A dataframe of results sorted by p-value, also containing effective sample size (ESS) and correlation estimate.
#' @importFrom SpatialPack modified.ttest
#' @seealso \link[SpatialPack]{modified.ttest}
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
wrapModTtest = function(X, Y, Cx, Ey, mapToFinest = FALSE, jointCoordinates = FALSE){
    n = nrow(X);m = nrow(Y)
    if(!jointCoordinates){
        matchedCoords = matchCoords(Cx, Ey, mapToFinest = mapToFinest)
        coordMat = matchedCoords$coordMat
        if(matchedCoords$mapToX){
            X = X[matchedCoords$id,]
        } else {
            Y = Y[matchedCoords$id,]
        }
    } else {
        coordMat = Cx
    }
    featGrid = expand.grid("featX" = colnames(X), "featY" = colnames(Y))
    out = simplify2array(loadBalanceBplapply(seq_len(nrow(featGrid)), function(i){
        unlist(modified.ttest(X[, featGrid[i, "featX"]], Y[, featGrid[i, "featY"]],
                       coordMat)[c("corr", "ESS", "p.value")])
    }))
    colnames(out) = makeNames(colnames(X), colnames(Y))
    rownames(out) = c("Correlation", "Effective sample size", "pVal")
    t(out)
}
#' Match coordinates to the nearest neighbour
#'
#' @inheritParams sbivarSingle
#'
#' @returns A list with components
#' \item{coordMat}{The new, single coordinate matrix}
#' \item{matToX}{A boolean, should X coordinates be taken}
#' \item{id}{The index of the coordinates to retain}
#' @importFrom RANN nn2
matchCoords = function(Cx, Ey, mapToFinest){
    if(mapToX <- xor(nrow(Cx) > nrow(Ey), mapToFinest)){
        id <- nn2(Cx, Ey, k = 1)$nn.idx[, 1]
        coordMat = Cx[id,]
    } else {
        id <- nn2(Ey, Cx, k = 1)$nn.idx[, 1]
        coordMat = Ey[id,]
    }
    list("coordMat" = coordMat, "mapToX" = mapToX, "id" = id)
}
