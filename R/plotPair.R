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
#' @param x,y Outcome vectors
#'
#' @returns A ggplot2 object
#' @export
#' @examples
#' example(sbivarSingle, "sbivar")
#' plotPairSingle(X[,1], Y[, 1], Cx, Ey, features = c("X1", "Y1"))
#' plotTopResultsSingle(resModtTest, Cx, Ey)
#' plotPairMulti(Xl, Yl, Cxl, Eyl, featx = "X1", featy = "Y1")
#' plotTopResultsMulti(estCorrelations, Xl, Yl, Cxl, Eyl)
plotPairSingleVectors = function(x, y, Cx, Ey, modalityNames = c("feature 1", "feature 2"), size = 2){
    stopifnot(length(x)==nrow(Cx), length(y)==nrow(Ey), ncol(Ey)==2, ncol(Cx)==2)
    coordMat = rbind(Cx, Ey);colnames(coordMat) = c("x", "y")
    plotDf = data.frame("outcome" = c(scaleZeroOne(x), scaleZeroOne(y)), coordMat,
                        "feature" = rep(modalityNames, times = c(length(x), lenght(y))))
    ggplot(data = plotDf, aes(x=x, y=y, col = outcome)) +
        geom_point(size = size) + facet_grid(feature~.) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("Dimension 1") + ylab("Dimension 2") + coord_fixed() +
        theme(axis.text = element_blank())
}
#' @rdname plotPairSingle
plotTopResultsSingle = function(resultsSingle, X, Y, Cx, Ey, normFun = logNorm,
                                ...){
    foo = checkInputSingle(X, Y, Cx, Ey)
    topPair = rownames(resultsSingle)[1]
    topFeats = strsplit(topPair, split = "--")[[1]]
    plotPairSingleVectors(x = normFun(X)[, topFeats[1]], y = logNorm(Y)[, topFeats[2]],
                   Cx = Cx, Ey = Ey)

}
plotPairSingle = function(X, Y, Cx, Ey, normFun = logNorm, features, ...){
    foo = checkInputSingle(X, Y, Cx, Ey)
    plotPairSingleVectors(x = normFun(X)[, features[1]], y = logNorm(Y)[, features[2]],
                   Cx = Cx, Ey = Ey)

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
