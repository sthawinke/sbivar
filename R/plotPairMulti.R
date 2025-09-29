#' Plot a chosen or top feature in the analysis for a series of images
#'
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
#' @param resultsMulti The results list, from call to \link{extractResultsMulti}
#' @inheritParams sbivarMulti
#' @inheritParams plotPairSingle
#' @examples
#' example(fitLinModels, "sbivar")
#' # Plot the feature pair with the most significant signal
#' plotTopResultsMulti(resGams, Xl, Yl, Cxl, Eyl)
#' # Plot an arbitrary feature pair
#' plotPairMulti(Xl, Yl, Cxl, Eyl, features = c("X1", "Y1"))
#' @export
#' @seealso \link{extractResultsMulti}, \link{sbivarMulti}, \link{fitLinModels}
plotTopResultsMulti = function(resultsMulti, Xl, Yl, Cx, Ey, parameter = "Intercept"){
    stopifnot(parameter %in% names(resultsMulti))
    topPair = rownames(resultsMulti[[parameter]])[1]
    topFeats = sund(topPair)
    plotPairMulti(Xl, Yl, Cxl, Eyl, features = topFeats)
}
#' @rdname plotTopResultsMulti
#' @export
#' @param features Feature vector of length 2 to be plotted
plotPairMulti = function(Xl, Yl, Cxl, Eyl, features, normalizationX = c("none", "log"),
                         normalizationY = c("none", "log"), modalityNames = c("Modality X", "Modality Y"),
                         size = 2){
    foo = checkInputMulti(Xl, Yl, Cxl, Eyl)
    stopifnot(length(features)==2)
    normalizationX = match.arg(normalizationX);normalizationY = match.arg(normalizationY)
    normFunX = switch(normalizationX, "none" = identity, "log" = logNorm)
    normFunY = switch(normalizationY, "none" = identity, "log" = logNorm)
    theme_set(theme_bw())
    dfList = Reduce(f = rbind, lapply(names(Xl), function(nam){
        coordMat = rbind(Cxl[[nam]], Eyl[[nam]]);colnames(coordMat) = c("x", "y")
        data.frame("outcome" = c(scaleHelpFun(Xl[[nam]], feat = features[1], normFun = normFunX),
                                 scaleHelpFun(Yl[[nam]], feat = features[2], normFun = normFunY)),
                   "image" = nam, coordMat,
                   "feature" = rep(modalityNames, times = c(nrow(Xl[[nam]]), nrow(Yl[[nam]]))))
    }))
    ggplot(data = dfList, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size) + facet_grid(image~feature) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("Dimension 1") + ylab("Dimension 2") + coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}

