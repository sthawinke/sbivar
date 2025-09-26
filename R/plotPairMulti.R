#' Plot a chosen or top feature in the analysis for a series of images
#'
#' @inheritParams sbivarMulti
#' @examples
#' plotPairMulti(Xl, Yl, Cxl, Eyl, featx = "X1", featy = "Y1")
#' plotTopResultsMulti(estCorrelations, Xl, Yl, Cxl, Eyl)
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
