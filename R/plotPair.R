#' Plot a top feature pair from the analysis, or an arbitrary pair for a series of images
#' @description Plot a chosen feature pair, or the highest ranking feature pair,
#' for a single image or multiple images.
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
#' The default is the intercept, i.e. the overall effect.
#' @param results The results list, from call to \link{extractResultsMulti}
#' @param rank An integer, the feature pair with the rank-th smallest p-value is plotted
#' @param ... passed onto lower level functions
#' @inheritParams sbivarMulti
#' @inheritParams plotPairSingle
#' @seealso \link{extractResultsMulti}, \link{fitLinModels}
#' @return A ggplot object
#' @order 1
#' @export
#' @details For sequence count data, such as transcriptomics, log-normalization
#' may be indicated to achieve clear plots (normalizationX = "log", see  \link{logNorm}). The normalization used for plotting is
#'  not necessarily the same as the one used for the analysis.
#' @examples
#' ### Single image
#' example(sbivarSingle, "sbivar")
#' # Plot the feature pair with the most significant signal
#' plotTopPair(resModtTest, X, Y, Cx, Ey)
#' # Plot an arbitrary feature pair
#' plotPairSingle(X, Y, Cx, Ey, features = c("X1", "Y1"))
#' ### Multi image
#' example(fitLinModels, "sbivar")
#' # Plot the feature pair with the most significant signal for a certain parameter,
#' #here "cofactor"
#' plotTopPair(resGams, Xl, Yl, Cxl, Eyl, parameter = "cofactor")
#' # Plot an arbitrary feature pair
#' plotPairMulti(Xl, Yl, Cxl, Eyl, features = c("X1", "Y1"))
plotTopPair = function(results, rank = 1, parameter = "Intercept", ...){
    stopifnot(is.numeric(rank))
    if(results$multiplicity == "single"){
        topFeats = sund(rownames(results$result)[rank])
        plotPairSingle(X = X, Y = Y, features = topFeats, Cx = Cx, Ey = Ey,
                       assayX = results$assayX, assayY = results$assayY, ...)
    } else if(results$multiplicity == "multi"){
        stopifnot(parameter %in% names(results$result))
        topFeats = sund(rownames(results$result)[[parameter]][rank])
plotPairMulti(Xl, Yl, Cxl, Eyl, features = topFeats,
              assayX = results$assayX, assayY = results$assayY, ...)
    } else {
        stop("Multiplicity tag unknown")
    }
}
#' @rdname plotTopPair
#' @export
#' @inheritParams plotPairSingle
#' @order 3
plotPairMulti = function(Xl, Yl, Cxl, Eyl, features, normalizationX = c("none", "log"),
                         normalizationY = c("none", "log"), size = 1.25, assayX, assayY){
    if(inherits(Xl[[1]], "SpatialExperiment")){
        Cxl = lapply(Xl, spatialCoords)
        Xl = lapply(Xl, assayT, assayX)
    }
    if(inherits(Yl[[1]], "SpatialExperiment")){
        Exl = lapply(Yl, spatialCoords)
        Yl = lapply(Yl, assayT, assayX)
    }
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
                   "feature" = rep(features, times = c(nrow(Xl[[nam]]), nrow(Yl[[nam]]))))
    }))
    ggplot(data = dfList, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size) + facet_grid(image~feature) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("Dimension 1") + ylab("Dimension 2") + coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}
#' @inheritParams sbivarSingle
#' @param results Results returned by \link{sbivarSingle}
#' @param x,y Outcome vectors
#' @param normalizationX,normalizationY Character strings, indicating what normalization is required
#' for X and Y matrices, respectively, before plotting, see details.
#' @param size Point size
#' @param features Feature vector of length 2 to be plotted
#' @rdname plotTopPair
#' @export
#' @order 2
plotPairSingle = function(X, Y, Cx, Ey, features, normalizationX = c("none", "log"),
                          normalizationY = c("none", "log"), ...){
    stopifnot(length(features)==2)
    if(inherits(X, "SpatialExperiment")){
        Cx = spatialCoords(Xl)
        X = assayT(X, assayX)
    }
    if(inherits(Y, "SpatialExperiment")){
        Ey = spatialCoords(Y)
        Y = assayT(Y, assayY)
    }
    foo = checkInputSingle(X, Y, Cx, Ey)
    normalizationX = match.arg(normalizationX);normalizationY = match.arg(normalizationY)
    normFunX = switch(normalizationX, "none" = identity, "log" = logNorm)
    normFunY = switch(normalizationY, "none" = identity, "log" = logNorm)
    plotPairSingleVectors(x = scaleHelpFun(feat = features[1], normFun = normFunX, X = X),
                          y = scaleHelpFun(feat = features[2], normFun = normFunY, X = Y),
                          Cx = Cx, Ey = Ey, modalityNames = features, ...)
}
#' @rdname plotTopPair
#' @param modalityNames Names to be given to the modalities,
#' appearing in the strip text of the columns. For plotTopPair() and
#' plotPairSingle(), the feature names are used.
#' @order 4
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

