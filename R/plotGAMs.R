#' Plot the fitted splines, and the correlation between them
#'
#' Spawns a three-panel plot with the splines fitted for the two variables, plus a visualization
#' of regions with positive and negative correlations between those splines.
#' The splines are refitted using \link{fitGAM}, so no gam-objects can be provided.
#' plotGAMsFromMatrix() takes entire outcome matrices X and Y as argument,
#' to be able to account for offsets in refitting the GAMs.#'
#' plotGAMsTopResults() plots the feature with the smallest p-value in resultsSingle.
#'
#' @note Both spline surfaces are scaled to the [-1,1] range, the same as
#' the correlation has naturally, for legibility.
#'
#' @inheritParams wrapGAMs
#' @param x,y outcome vectors
#' @param newGrid The grid in which to evaluate the GAMs
#' @param offsets List of length two with offsets
#' @param scaleFun The scaling function to be applied before plotting
#' @param addTitle A boolean, should a title be plotted
#' @param features The features to plot
#' @param resultsSingle Result of a call to \link{sbivarSingle}
#' @param ... passed onto \link{fitGAM}
#'
#' @returns A ggplot object
#' @export
#'
#' @examples
#' example(sbivarSingle, "sbivar")
#' plotGAMs(X[, 1], Y[, 1], Cx, Ey)
#' plotGAMsFromMatrix(X, Y, features = c("X1", "Y1"), Cx = Cx, Ey = Ey)
#' plotGAMsTopResults(resGAMs, X, Y, Cx = Cx, Ey = Ey)
#' @import ggplot2
#' @importFrom reshape2 melt
#' @order 1
plotGAMs = function(x, y, Cx, Ey, newGrid, offsets = list(), scaleFun = "scaleMinusOne",
                    families = list("X" = gaussian(), "Y" = gaussian()),
                    addTitle = TRUE, n_points_grid = min(length(x), length(y)),
                    features = c("x", "y"), ...){
    stopifnot(is.numeric(n_points_grid), all(vapply(families, FUN.VALUE = TRUE, is, "family")))
    colnames(Cx) = colnames(Ey) = c("x", "y")
    scaleFun = get(as.character(scaleFun), mode = "function", getNamespace("sbivar"))
    if(missing(newGrid))
        newGrid = buildNewGrid(Cx = Cx, Ey = Ey, n_points_grid = n_points_grid)
    modelx <- fitGAM(df = data.frame("value" = x, Cx), outcome = "value",
                     family = families[["X"]], offset = offsets[["X"]], ...)
    modely <- fitGAM(df = data.frame("value" = y, Ey), family = families[["Y"]],
                     offset = offsets[["Y"]], outcome = "value", ...)
    predx = vcovPredGam(modelx, newdata = newGrid)
    predy = vcovPredGam(modely, newdata = newGrid)
    corContr = (predx$pred-mean(predx$pred))*(predy$pred-mean(predy$pred))
    corEst = sum(corContr)/((nrow(newGrid)-1)*sd(predx$pred)*sd(predy$pred))
    dat = rbind(data.frame(newGrid, value = scaleFun(predx$pred), feature = "x"),
                data.frame(newGrid, value = scaleFun(predy$pred), feature = "y"),
                data.frame(newGrid, value = scaleFun(corContr), feature = "cor"))
    gridMolt = melt(dat, id.vars = c("x", "y", "feature"), value.name = "Value")
    gridMolt$feature = factor(gridMolt$feature, levels = c("x", "y", "cor"),
                              labels = c(features[1], features[2], "Correlation"),
                              ordered = TRUE)
    ggplot(gridMolt, aes(x, y, fill = Value)) +
        geom_raster() + coord_fixed() +
        facet_grid( ~ feature) +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_fill_viridis_c(option = "H", name = "") +
        if(addTitle)
            labs(title = paste("Spline surfaces, and contributions to correlation estimate", round(corEst, 3)))
}
#' @export
#' @rdname plotGAMs
#' @order 3
plotGAMsFromMatrix = function(X, Y, features, Cx, Ey,
                    families = list("X" = gaussian(), "Y" = gaussian()), ...){
    if(families[["X"]]$family!="gaussian"){
        X = X[idX <- (rowSums(X)>0),]
        Cx = Cx[idX,]
    }
    if(families[["Y"]]$family!="gaussian"){
        Y = Y[idY <- (rowSums(Y)>0),]
        Ey = Ey[idY,]
    }
    offsets = list("X" = makeOffset(X, families[["X"]]),
                   "Y" = makeOffset(Y, families[["Y"]]))
    plotGAMs(x = X[, features[1]], y = Y[,features[2]], Cx = Cx, Ey = Ey,
             families = families, offsets = offsets, features = features, ...)
}
#' @export
#' @rdname plotGAMs
#' @order 2
#' @inheritParams plotTopPair
plotGAMsTopResults = function(resultsSingle, X, Y, Cx, Ey, topRank = 1, ...){
    topFeats = sund(rownames(resultsSingle$result)[topRank])
    plotGAMsFromMatrix(X = X, Y = Y, features = topFeats, Cx = Cx, Ey = Ey, ...)
}
#' Make a list of offsets
#'
#' @param X The data matrix
#' @param family The distribution family
#' @returns A list of length two with offsets
makeOffset = function(X, family){
    out = if(family$family != "gaussian"){
        libSizes = rowSums(X)
        switch(family$link, "inverse" = 1/libSizes, "log" = log(libSizes))
    }
    return(out)
}
