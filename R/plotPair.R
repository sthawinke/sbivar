#' Plot a top feature pair from the analysis, or an arbitrary pair for a series of images
#' @description Plot a chosen feature pair, or the highest ranking feature pair,
#' for a single image or multiple images.
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
#' The default is the intercept, i.e. the overall effect.
#' @param results The results list, from call to \link{extractResultsMulti}
#' @param topRank An integer, the feature pair with the rank-th smallest p-value is plotted
#' @param ... passed onto lower level functions
#' @inheritParams sbivarMulti
#' @inheritParams plotPairSingle
#' @seealso \link{extractResultsMulti}, \link{fitLinModels}
#' @return A ggplot object
#' @order 1
#' @export
#' @details For sequence count data, such as transcriptomics, log-normalization
#' may be indicated to achieve clear plots (normX = "log", see  \link{logNorm}). The normalization used for plotting is
#'  not necessarily the same as the one used for the analysis.
#' @examples
#' ### Single image
#' example(sbivar, "sbivar")
#' # Plot the feature pair with the most significant signal
#' plotTopPair(resMoransI, X, Y, Cx, Ey)
#' # Plot an arbitrary feature pair
#' plotPairSingle(X, Y, Cx, Ey, features = c("X1", "Y1"))
#' ### Multi image
#' example(fitLinModels, "sbivar")
#' # Plot the feature pair with the most significant signal for a certain parameter,
#' #here the intercept (overall effect)
#' plotTopPair(resMoran, Vicari$TranscriptOutcomes, Vicari$MetaboliteOutcomes,
#' Vicari$TranscriptCoords, Vicari$MetaboliteCoords, parameter = "Intercept")
#' # Plot an arbitrary feature pair
#' plotPairMulti(Vicari$TranscriptOutcomes, Vicari$MetaboliteOutcomes, normX = "log", normY = "log",
#' Vicari$TranscriptCoords, Vicari$MetaboliteCoords, features = c("mt.Nd2", "X555.20713"))
plotTopPair = function(results, ..., normX = results$normX, normY = results$normY, topRank = 1, parameter = "Intercept"){
    stopifnot(is.numeric(topRank))
    if(!results$multi){
        topFeats = sund(rownames(results$result)[topRank])
        plotPairSingle(features = topFeats, assayX = results$assayX,
                       assayY = results$assayY, normX = normX,
                       normY = normY,...)
    } else {
        stopifnot(parameter %in% names(results$result))
        topFeats = sund(rownames(results$result[[parameter]])[topRank])
        plotPairMulti(features = topFeats, assayX = results$assayX,
              assayY = results$assayY, normX = normX,
              normY = normY, ...)
    }
}
#' @rdname plotTopPair
#' @export
#' @inheritParams plotPairSingle
#' @order 3
#' @param theme the ggplot2 theme
plotPairMulti = function(Xl, Yl, Cxl, Eyl, features, normX = c("none", "log"),
                         normY = c("none", "log"), size = 1.25, assayX, assayY, theme = theme_bw()){
    Xl =  lapply(getX(Xl, assayX), addDimNames, "X")
    Yl = lapply(getX(Yl, assayY), addDimNames, "Y")
    Cxl = getSpatialCoords(Xl, Cxl)
    Eyl = getSpatialCoords(Yl, Eyl)
    foo = checkInputMulti(Xl, Yl, Cxl, Eyl)
    features = make.names(features)
    stopifnot(length(features)==2)
    normX = match.arg(normX);normY = match.arg(normY)
    theme_set(theme)
    nfx = getNormFun(normX);nfy = getNormFun(normY)
    dfList = Reduce(f = rbind, lapply(names(Xl), function(nam){
        X = nfx(Xl[[nam]]);Y= nfy(Yl[[nam]])
        coordMat = rbind(Cxl[[nam]][rownames(X),], Eyl[[nam]][rownames(Y),]);colnames(coordMat) = c("x", "y")
        data.frame("outcome" = c(scaleHelpFun(X, feat = features[1]),
                                 scaleHelpFun(Y, feat = features[2])),
                   "image" = nam, coordMat,
                   "feature" = rep(features, times = c(nrow(Xl[[nam]]), nrow(Yl[[nam]]))))
    }))
    ggplot(data = dfList, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size) + facet_grid(image~feature) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("x coordinate") + ylab("y coordinate") + coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}
#' @inheritParams sbivar
#' @param results Results returned by \link{sbivarSingle}
#' @param x,y Outcome vectors
#' @param normX,normY Character strings, indicating what normalization is required
#' for X and Y matrices, respectively, before plotting, see details.
#' @param size Point size
#' @param features Feature vector of length 2 to be plotted
#' @rdname plotTopPair
#' @export
#' @order 2
plotPairSingle = function(X, Y, Cx, Ey, features, normX = c("none", "log"),
                          normY = c("none", "log"), assayX, assayY, ...){
    stopifnot(length(features)==2)
    if(inherits(X, "SpatialExperiment")){
        Cx = spatialCoords(X)
        X = assayT(X, assayX)
    }
    if(inherits(Y, "SpatialExperiment")){
        Ey = spatialCoords(Y)
        Y = assayT(Y, assayY)
    }
    X = addDimNames(X, "X");Y = addDimNames(Y, "Y")
    features = make.names(features)
    foo = checkInputSingle(X, Y, Cx, Ey)
    normX = match.arg(normX);normY = match.arg(normY)
    rownames(Cx) = rownames(X);rownames(Ey) = rownames(Y)
    X = getNormFun(normX)(X)
    Y = getNormFun(normY)(Y)
    plotPairSingleVectors(x = scaleHelpFun(feat = features[1], X),
                          y = scaleHelpFun(feat = features[2], Y),
                          Cx = Cx[rownames(X),], Ey = Ey[rownames(Y),], modalityNames = features, ...)
}
#' @rdname plotTopPair
#' @param modalityNames Names to be given to the modalities,
#' appearing in the strip text of the columns. For plotTopPair() and
#' plotPairSingle(), the feature names are used.
#' @order 4
plotPairSingleVectors = function(x, y, Cx, Ey, size = 1.25,
                                 modalityNames = c("Modality X", "Modality Y"), theme = theme_bw(),... ){
    theme_set(theme)
    stopifnot(length(x)==nrow(Cx), length(y)==nrow(Ey), ncol(Ey)==2, ncol(Cx)==2)
    coordMat = rbind(Cx, Ey);colnames(coordMat) = c("x", "y")
    plotDf = data.frame("outcome" = c(x, y), coordMat,
                        "feature" = rep(modalityNames, times = c(length(x), length(y))))
    ggplot(data = plotDf, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size) + facet_grid(~feature) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("x coordinate") + ylab("y coordinate") + coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}

