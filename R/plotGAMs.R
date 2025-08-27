#' Title
#'
#' @inheritParams wrapGAMs
#' @inheritParams testManyGAMs
#' @param offsets List of length two with offsets
#' @param scaleFun The scaling function to be applied before plotting
#' @param addTitle A booleam, should a title be plotted
#' @param ... passed onto \link{fitGAM}
#'
#' @returns A ggplot object
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom reshape2 melt
plotGAMs = function(x, y, Cx, Ey, newGrid, families, offsets,
                       scaleFun = "scaleMinusOne", addTitle = TRUE, ...){
    scaleFun = match.fun(scaleFun)
    if(missing(newGrid))
        newGrid = buildNewGrid(Cx = Cx, Ey = Ey, n_points_grid = n_points_grid)
    modelx <- fitGAM(df = data.frame("value" = x, Cx), family = families[["X"]],
                     offset = offsets[["X"]], ...)
    modely <- fitGAM(df = data.frame("value" = y, Ey), family = families[["Y"]],
                     offset = offsets[["Y"]], ...)
    predx = vcovPredGam(modelx, newdata = newGrid)
    predy = vcovPredGam(modely, newdata = newGrid)
    corContr = (cen1 <- (predx$pred-mean(predx$pred)))*(cen2 <- (predy$pred-mean(predy$pred)))
    corEst = sum(corContr)/((nrow(newGrid)-1)*sd(predx$pred)*sd(predy$pred))
    dat = rbind(data.frame(newGrid, value = scaleFun(predx$pred), feature = "x"),
                data.frame(newGrid, value = scaleFun(predy$pred), feature = "y"),
                data.frame(newGrid, value = scaleFun(corContr), feature = "cor"))
    gridMolt = melt(dat, id.vars = c("x", "y", "feature"), value.name = "Value")
    gridMolt$feature = factor(gridMolt$feature, levels = c("x", "y", "cor"),
                              labels = c("x", "y", "cor"), ordered = TRUE)
    ggplot(gridMolt, aes(x, y, fill = Value)) +
        geom_raster() + coord_fixed() +
        facet_grid( ~ feature) +
        scale_fill_viridis_c(option = "H", name = "") +
        if(addTitle)
            labs(title = paste("Estimated spline surfaces, and contributions to correlation estimate", round(corEst, 3)))
}
#' A wrapper frunction to plot from matrices
#'
#' @inheritParams sbivar
#' @param featx,featy The features to plot
#' @param ... Passed onto \link{plotGAMs}.
#'
#' @returns See \link{plotGAMs}
#' @export
#'
#' @examples
plotGAMsFromMatrix = function(X, Y, featx, featy, Cx, Ey, families, ...){
    offsets = list("X" = makeOffset(X, families[["X"]]),
                   "Y" = makeOffset(Y, families[["Y"]]))
    plotGAMs(x = X[, featx], Y = Y[,featY], Cx = Cx, Ey = Ey,
             families = families, offsets = offsets, ...)
}
#' Make a list of offsets
#'
#' @param X The data matrix
#' @param family The distribution family
#'
#' @returns A list of length two with offsets
makeOffset = function(X, family){
    out = if(family$family != "gaussian"){
        libSizes = rowSums(X)
        switch(family$link, "inverse" = 1/libSizes, "log" = log(libSizes))
    }
    return(out)
}
