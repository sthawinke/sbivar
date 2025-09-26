#' Plot a feature pair
#'
#' Plot a chosen feature pair, or the highest ranking feature pair, for a singlepair of images,
#'
#' @details For sequence count data, such as transcriptomics, log-normalization
#' seems indicated to achieve clear plots.
#'
#' @inheritParams sbivarSingle
#' @param resultsSingle,resultsMulti Results, returned by \link{sbivarSingle} and
#' \link{extractResultsMulti}, respectively.
#' @param featx,featy The features in X and Y to be plotted
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
#' @param x,y Outcome vectors
#' @param normalizationX,normalizationY A character string, indicating what normalization is required
#' for X and Y matrrices, respectively, before plotting.
#' This is not the same normalization as used for the analysis, see details.
#' @param size Point size
#' @param ... Passed onto \link{plotPairSingleVectors}
#'
#' @returns A ggplot2 object
#' @export
#' @examples
#' example(sbivarSingle, "sbivar")
#' plotPairSingle(X, Y, Cx, Ey, features = c("X1", "Y1"))
#' plotTopResultsSingle(resModtTest, X, Y, Cx, Ey)
plotTopResultsSingle = function(resultsSingle, X, Y, Cx, Ey, ...){
    topPair = rownames(resultsSingle)[1]
    topFeats = strsplit(topPair, split = "__")[[1]]
    plotPairSingle(X = X, Y = Y, features = topFeats, Cx = Cx, Ey = Ey, ...)

}
plotPairSingle = function(X, Y, Cx, Ey, features, normalizationX = c("none", "log"),
                          normalizationY = c("none", "log"), ...){
    normalizationX = match.arg(normalizationX);normalizationY = match.arg(normalizationY)
    normFunX = switch(normalizationX, "none" = identity, "log" = logNorm)
    normFunY = switch(normalizationY, "none" = identity, "log" = logNorm)
    foo = checkInputSingle(X, Y, Cx, Ey)
    plotPairSingleVectors(x = normFunX(X)[, features[1]], y = normFunY(Y)[, features[2]],
                   Cx = Cx, Ey = Ey, ...)

}
#' @rdname plotTopResultsSingle
plotPairSingleVectors = function(x, y, Cx, Ey, size = 2,
                                 modalityNames = c("feature 1", "feature 2")){
    theme_set(theme_bw())
    stopifnot(length(x)==nrow(Cx), length(y)==nrow(Ey), ncol(Ey)==2, ncol(Cx)==2)
    coordMat = rbind(Cx, Ey);colnames(coordMat) = c("x", "y")
    plotDf = data.frame("outcome" = c(scaleZeroOne(x), scaleZeroOne(y)), coordMat,
                        "feature" = rep(modalityNames, times = c(length(x), length(y))))
    ggplot(data = plotDf, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size) + facet_grid(~feature) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("Dimension 1") + ylab("Dimension 2") + coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}
