#' Plot a specified feature pair for the data provided
#' @rdname plotPair
#' @export
#' @inheritParams plotPairSingle
#' @inheritParams sbivarMulti
#' @order 2
#' @param theme the ggplot2 theme
plotPairMulti <- function(Xl, Yl, Cxl, Eyl, features, normX = c("none", "rel", "log"), scaleBySampleSums = FALSE,
    normY = c("none", "rel", "log"), size = 1.25, assayX, assayY, theme = theme_bw()) {
    stopifnot(is.logical(scaleBySampleSums), is.numeric(size))
    Xl <- getX(Xl, assayX)
    Yl <- getX(Yl, assayY)
    Cxl <- getSpatialCoords(Xl, Cxl)
    Eyl <- getSpatialCoords(Yl, Eyl)
    foo <- checkInputMulti(Xl, Yl, Cxl, Eyl)
    features <- make.names(features)
    stopifnot(length(features) == 2)
    normX <- match.arg(normX)
    normY <- match.arg(normY)
    theme_set(theme)
    dfList <- do.call(rbind, lapply(names(Xl), function(nam) {
        size <- getSize(Xl[[nam]], Yl[[nam]], normX, normY, scaleBySampleSums, size = size)
        X <- normMat(Xl[[nam]], normX)
        Y <- normMat(Yl[[nam]], normY)
        coordMat <- rbind(Cxl[[nam]][rownames(X), ], Eyl[[nam]][rownames(Y), ])
        colnames(coordMat) <- c("x", "y")
        data.frame(
            "outcome" = c(
                scaleHelpFun(X, feat = features[1]),
                scaleHelpFun(Y, feat = features[2])
            ),
            "size" = size[c(rownames(X), rownames(Y))], "image" = nam, coordMat,
            "feature" = rep(features, times = c(nrow(X), nrow(Y)))
        )
    }))
    ggplot(data = dfList, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size) +
        facet_grid(image ~ feature) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("x coordinate") +
        ylab("y coordinate") +
        coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}
#' @inheritParams sbivar
#' @param x,y Outcome vectors
#' @param normX,normY Character strings, indicating what normalization is required
#' for X and Y matrices, respectively, before plotting, see details.
#' @param size Point size
#' @param features Feature vector of length 2 to be plotted
#' @param scaleBySampleSums boolean, should the point size be scaled by sample sums?
#' @rdname plotPair
#' @export
#' @order 1
plotPairSingle <- function(
      X, Y, Cx, Ey, features, normX = c("none", "rel", "log"),
      normY = c("none", "rel", "log"), assayX, assayY, scaleBySampleSums = FALSE, size = 1.5, ...
) {
    stopifnot(length(features) == 2, is.numeric(size), is.logical(scaleBySampleSums))
    if (inherits(X, "SpatialExperiment")) {
        Cx <- SpatialExperiment::spatialCoords(X)
        X <- assayT(X, assayX)
    }
    if (inherits(Y, "SpatialExperiment")) {
        Ey <- SpatialExperiment::spatialCoords(Y)
        Y <- assayT(Y, assayY)
    }
    features <- make.names(features)
    foo <- checkInputSingle(X, Y, Cx, Ey)
    normX <- match.arg(normX)
    normY <- match.arg(normY)
    size <- getSize(X, Y, normX, normY, scaleBySampleSums, size = size)
    X <- normMat(X, normX)
    Y <- normMat(Y, normY)
    plotPairSingleVectors(
        x = scaleHelpFun(feat = features[1], X), size = size[c(rownames(X), rownames(Y))],
        y = scaleHelpFun(feat = features[2], Y),
        Cx = Cx[rownames(X), ], Ey = Ey[rownames(Y), ], modalityNames = features, ...
    )
}
#' @rdname plotPair
#' @param modalityNames Names to be given to the modalities,
#' appearing in the strip text of the columns. For plotTopPair() and
#' plotPairSingle(), the feature names are used.
#' @order 3
plotPairSingleVectors <- function(x, y, Cx, Ey, size,
    modalityNames = c("Modality X", "Modality Y"), theme = theme_bw(), ...) {
    theme_set(theme)
    stopifnot(length(x) == nrow(Cx), length(y) == nrow(Ey), ncol(Ey) == 2, ncol(Cx) == 2)
    coordMat <- rbind(Cx, Ey)
    colnames(coordMat) <- c("x", "y")
    plotDf <- data.frame(
        "outcome" = c(x, y), coordMat,
        "feature" = factor(rep(modalityNames, times = c(length(x), length(y))),
            levels = modalityNames, ordered = TRUE
        ),
        "size" = size
    )
    ggplot(data = plotDf, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size) +
        facet_grid(~feature) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("x coordinate") +
        ylab("y coordinate") +
        coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}
