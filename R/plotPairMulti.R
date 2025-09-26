#' Plot a chosen or top feature in the analysis for a series of images
#'
#' @inheritParams sbivarMulti
#' @examples
#' example(fitLinModels, "sbivar")
#' plotPairMulti(Xl, Yl, Cxl, Eyl, featx = "X1", featy = "Y1")
#' plotTopResultsMulti(resGams, Xl, Yl, Cxl, Eyl)
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
plotTopResultsMulti = function(resultsMulti, Xl, Yl, Cx, Ey, parameter = "Intercept"){
    stopifnot(parameter %in% names(resultsMulti))
    topPair = rownames(resultsMulti[[parameter]][1])
    topFeats = sund(topPair)
    plotPairMulti(Xl, Yl, Cxl, Eyl, features = topFeats)
}
#' @rdname plotPairSingle

plotPairMulti = function(Xl, Yl, Cxl, Eyl, features, normalizationX = c("none", "log"),
                         normalizationY = c("none", "log")){
    normalizationX = match.arg(normalizationX);normalizationY = match.arg(normalizationY)
    normFunX = switch(normalizationX, "none" = identity, "log" = logNorm)
    normFunY = switch(normalizationY, "none" = identity, "log" = logNorm)
    foo = checkInputMulti(Xl, Yl, Cxl, Eyl)
    theme_set(theme_bw())
    dfList = mapply(Xl, Yl, Cxl, Eyl, FUN = function(X, Y, Cx, Ey){
        coordMat = rbind(Cx, Ey);colnames(coordMat) = c("x", "y")
        plotDf = data.frame("outcome" = c(scaleZeroOne(x), scaleZeroOne(y)), coordMat,
                            "feature" = rep(modalityNames, times = c(length(x), length(y))))
    })

    ggplot(data = plotDf, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size) + facet_grid(~feature) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("Dimension 1") + ylab("Dimension 2") + coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}
makePlotSingleDf = function(gene, met, sample, normFun){
    magpieCoordsShifted = lapply(magpieCoords, shiftCoord)
    locObjsStxWarped = lapply(magpieCoordsShifted[sample], function(x) x$Stx)
    locObjsMsiWarped = lapply(magpieCoordsShifted[sample], function(x) x$Msi)
    geneLocs = Reduce(locObjsStxWarped, f = rbind)
    metLocs = Reduce(locObjsMsiWarped, f = rbind)
    geneSn = rep(sample, times = sapply(locObjsStxWarped, nrow))
    metSn = rep(sample, times = sapply(locObjsMsiWarped, nrow))
    geneVal = unlist(lapply(exprObjs[sample], function(xx) if(gene %in% colnames(xx)) normFun(xx)[, gene] else rep(NA, nrow(xx))))
    metVal = unlist(lapply(msiObjs[sample], function(xx) if(met %in% colnames(xx)) normFun(xx)[, met] else rep(NA, nrow(xx))))
    plotDf = data.frame(rbind(data.frame(geneLocs, "image" = geneSn, feat = "Gene"),
                              data.frame(metLocs, "image" = metSn, feat = "Metabolite")),
                        "value" = c(scaleZeroOne(geneVal), scaleZeroOne(metVal)))
    plotDf$feat = factor(plotDf$feat, labels = c("Gene", "Metabolite"),
                         levels = c("Gene", "Metabolite"), ordered = TRUE)
    plotDf
}
