#' Plot a feature pair
#'
#' Plot a chosen feature pair, or the highest ranking feature pair, for a single pair of images
#'
#' @details For sequence count data, such as transcriptomics, log-normalization
#' seems indicated to achieve clear plots. The normalization used for plotting is
#'  not necessarily the same as the one used for the analysis.
#'
#' @inheritParams sbivarSingle
#' @param resultsSingle Results returned by \link{sbivarSingle}
#' @param x,y Outcome vectors
#' @param normalizationX,normalizationY Character strings, indicating what normalization is required
#' for X and Y matrices, respectively, before plotting, see details.
#' @param size Point size
#' @param ... Passed onto \link{plotPairSingleVectors}
#'
#' @returns A ggplot2 object
#' @export
#' @examples
#' example(sbivarSingle, "sbivar")
#' # Plot the feature pair with the most significant signal
#' plotTopResultsSingle(resModtTest, X, Y, Cx, Ey)
#' # Plot an arbitrary feature pair
#' plotPairSingle(X, Y, Cx, Ey, features = c("X1", "Y1"))
plotTopResultsSingle = function(resultsSingle, X, Y, Cx, Ey, ...){
    topPair = rownames(resultsSingle)[1]
    topFeats = sund(topPair)
    plotPairSingle(X = X, Y = Y, features = topFeats, Cx = Cx, Ey = Ey, ...)

}
#' @rdname plotTopResultsSingle
#' @export
#' @param features Feature vector of length 2 to be plotted
plotPairSingle = function(X, Y, Cx, Ey, features, normalizationX = c("none", "log"),
                          normalizationY = c("none", "log"), ...){
    stopifnot(length(features)==2)
    foo = checkInputSingle(X, Y, Cx, Ey)
    normalizationX = match.arg(normalizationX);normalizationY = match.arg(normalizationY)
    normFunX = switch(normalizationX, "none" = identity, "log" = logNorm)
    normFunY = switch(normalizationY, "none" = identity, "log" = logNorm)
    plotPairSingleVectors(x = scaleHelpFun(feat = features[1], normFun = normFunX, X = X),
                          y = scaleHelpFun(feat = features[2], normFun = normFunY, X = Y),
                   Cx = Cx, Ey = Ey, modalityNames = features,...)

}
#' @rdname plotTopResultsSingle
#' @param modalityNames Names to be given to the modalities,
#' appearing in the strip text of the columns
plotPairSingleVectors = function(x, y, Cx, Ey, size = 1.25,
                                 modalityNames = c("Modality X", "Modality Y")){
    theme_set(theme_bw())
    stopifnot(length(x)==nrow(Cx), length(y)==nrow(Ey), ncol(Ey)==2, ncol(Cx)==2)
    coordMat = rbind(Cx, Ey);colnames(coordMat) = c("x", "y")
    plotDf = data.frame("outcome" = c(x, y), coordMat,
                        "feature" = rep(modalityNames, times = c(length(x), length(y))))
    ggplot(data = plotDf, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size) + facet_grid(~feature) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("Dimension 1") + ylab("Dimension 2") + coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}
