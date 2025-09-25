#' Plot a feature pair
#'
#' Plot a chosen feature pair for a single or list of matrices,
#' or plot the highest ranking feature pair.
#'
#' @inheritParams sbivarSingle
#' @inheritParams sbivarMulti
#' @param resultsSingle,resultsMulti Results, returned by \link{sbivarSingle} and
#' \link{extractResultsNulti}, respectively.
#'
#' @returns A ggplot2 object
#' @export
#' @examples
plotPairSingle = function(x, y, Cx, Ey){

}
plotTopResultsSingle = function(resultsSingle, Cx, Ey){

}
plotPairMulti = function(Xl, Yl, Cxl, Eyl, featx, featy){
    foo = checkInputMultiple(Xl, Yl, Cxl, Eyl)
}
plotTopResultsMulti = function(resultsMulti, Xl, Yl, Cx, Ey,
                               parameter = "Intercept"){
    plotPairMulti(Xl, Yl, Cxl, Eyl, featx = featx, featy = featy)
}
