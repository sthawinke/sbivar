setMethod(
    "plotGAMsTopPair",
    "sbivarResultsGAMs",
    function(x, topRank, parameter, ...) {
        stopifnot(is.numeric(topRank))
        topFeats <- (
            if (x@multi) {
                x@result[[parameter]]
            } else {
                x@result
            })[topRank, c("Modality_X", "Modality_Y")]
        Cx <- getSpatialCoords(X, Cx)
        X <- getX(X, x@assayX)
        Ey <- getSpatialCoords(Y, Ey)
        Y <- getX(Y, x@assayY)
        plotGAMs(
            X = X, Y = Y, features = topFeats, Cx = Cx, Ey = Ey, families = x@families,
            multi = x@multi, normX = x@normX, normY = x@normY, Gamm = !x@multi && x@Gamm, correlation = x@correlation, ...
        )
    }
)
plotGAMsTopResults <- function(results, X, Y, Cx, Ey, topRank = 1,
                               parameter = "Intercept", families = results$families, ...) {

}
