#' Plot a feature pair
#'
#' Plot a chosen feature pair for a single or list of matrices,
#' or plot the highest ranking feature pair.
#'
#' @inheritParams sbivarSingle
#' @inheritParams sbivarMulti
#' @param resultsSingle,resultsMulti Results, returned by \link{sbivarSingle} and
#' \link{extractResultsMulti}, respectively.
#' @param featx,featy The features in X and Y to be plotted
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
#'
#' @returns A ggplot2 object
#' @export
#' @examples
#' example(sbivarSingle, sbivar)
#' plotPairSingle(X[,1], Y[, 1], Cx, Ey)
#' plotTopResultsSingle(resModtTest, Cx, Ey)
#' plotPairMulti(Xl, Yl, Cxl, Eyl, featx = "X1", featy = "Y1")
#' plotTopResultsMulti(estCorrelations, Xl, Yl, Cxl, Eyl)
plotPairSingle = function(x, y, Cx, Ey){

}
#' @rdname plotPairSingle
plotTopResultsSingle = function(resultsSingle, Cx, Ey){

}
#' @rdname plotPairSingle
plotPairMulti = function(Xl, Yl, Cxl, Eyl, featx, featy){
    foo = checkInputMulti(Xl, Yl, Cxl, Eyl)
}
#' @rdname plotPairSingle
plotTopResultsMulti = function(resultsMulti, Xl, Yl, Cx, Ey,
                               parameter = "Intercept"){
    plotPairMulti(Xl, Yl, Cxl, Eyl, featx = featx, featy = featy)
}
