#' Plot a chosen or top feature in the analysis for a series of images
#'
#' @inheritParams sbivarMulti
#' @inheritParams plotPairSingle
#' @examples
#' example(fitLinModels, "sbivar")
#' plotTopResultsMulti(resGams, Xl, Yl, Cxl, Eyl)
#' plotPairMulti(Xl, Yl, Cxl, Eyl, feauresx = c("X1", "Y1"))
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
plotTopResultsMulti = function(resultsMulti, Xl, Yl, Cx, Ey, parameter = "Intercept"){
    stopifnot(parameter %in% names(resultsMulti))
    topPair = rownames(resultsMulti[[parameter]][1])
    topFeats = sund(topPair)
    plotPairMulti(Xl, Yl, Cxl, Eyl, features = topFeats)
}
#' @rdname plotPairSingle
plotPairMulti = function(Xl, Yl, Cxl, Eyl, features, normalizationX = c("none", "log"),
                         normalizationY = c("none", "log"), size = 2){
    normalizationX = match.arg(normalizationX);normalizationY = match.arg(normalizationY)
    normFunX = switch(normalizationX, "none" = identity, "log" = logNorm)
    normFunY = switch(normalizationY, "none" = identity, "log" = logNorm)
    foo = checkInputMulti(Xl, Yl, Cxl, Eyl)
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

