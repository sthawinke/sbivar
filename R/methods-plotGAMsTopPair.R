#' @rdname plotGAMsTopPair
#' @examples
#' example(sbivar, "sbivar")
#' plotGAMsTopPair(resGAMs, X, Y, Cx, Ey)
#' @return a ggplot2 object
#' @param x a sbivarResultsGAMs object
setMethod(
    "plotGAMsTopPair",
    "sbivarResultsGAMs",
    function(x, X, Y, Cx, Ey, topRank, parameter, ...) {
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
