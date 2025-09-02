#' Plot the fitted splines, and the correlation between them
#'
#' @inheritParams wrapGAMs
#' @param x,y outcome vetors
#' @param offsets List of length two with offsets
#' @param scaleFun The scaling function to be applied before plotting
#' @param addTitle A booleam, should a title be plotted
#' @param ... passed onto \link{fitGAM}
#'
#' @returns A ggplot object
#' @export
#'
#' @examples
#' example(sbivarSingle, "sbivar")
#' plotGAMs(X[, 1], Y[, 1], Cx, Ey)
#' @import ggplot2
#' @importFrom reshape2 melt
plotGAMs = function(x, y, Cx, Ey, newGrid, offsets = list(), scaleFun = "scaleMinusOne",
                    families = list("X" = gaussian(), "Y" = gaussian()),
                    addTitle = TRUE, n_points_grid = 5e2, ...){
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
#' @inheritParams sbivarSingle
#' @param featX,featY The features to plot
#' @param ... Passed onto \link{plotGAMs}.
#'
#' @returns See \link{plotGAMs}
#' @export
#'
#' @examples
#' example(sbivarSingle, "sbivar")
#' plotGAMsFromMatrix(X, Y, featX = "X1", featY = "Y1", Cx = Cx, Ey = Ey)
plotGAMsFromMatrix = function(X, Y, featX, featY, Cx, Ey,
                              families = list("X" = gaussian(), "Y" = gaussian()), ...){
    offsets = list("X" = makeOffset(X, families[["X"]]),
                   "Y" = makeOffset(Y, families[["Y"]]))
    plotGAMs(x = X[, featX], y = Y[,featY], Cx = Cx, Ey = Ey,
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
