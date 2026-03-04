#' @title Plot a feature pair
#' @description Plot a chosen feature pair, or the highest ranking feature pair,
#' for a single image or multiple images.
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
#' The default is the intercept, i.e. the overall effect.
#' @param topRank An integer, the feature pair with the rank-th smallest p-value is plotted
#' @param ... passed onto lower level functions
#' @inheritParams sbivarMulti
#' @inheritParams plotPairSingle
#' @seealso \link{extractResultsMulti}, \link{fitLinModels}
#' @return A ggplot object
#' @order 1
#' @export
#' @details For sequence count data, such as transcriptomics, normalization
#' may be indicated to achieve clear plots (normX = "rel" or "log", see  \link{normMat}).
#' The normalization used for plotting is not necessarily the same as the one used for the analysis.
#' @examples
#' ### Single image
#' # Single image analysis on synthetic data
#' n <- 8e1
#' m <- 9e1
#' p <- 3
#' k <- 2
#' X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
#' Y <- matrix(rnorm(m * k), m, k, dimnames = list(NULL, paste0("Y", seq_len(k))))
#' Cx <- matrix(runif(n * 2), n, 2)
#' Ey <- matrix(runif(m * 2), m, 2)
#' colnames(Cx) <- colnames(Ey) <- c("x", "y")
#' resMoransI <- sbivar(X, Y, Cx, Ey, method = "Moran's I")
#' # Plot the feature pair with the most significant signal
#' plotTopPair(resMoransI, X, Y, Cx, Ey)
#' # Plot an arbitrary feature pair
#' plotPairSingle(X, Y, Cx, Ey, features = c("X1", "Y1"))
#' ### Multi image
#' data(Vicari)
#' # Plot an arbitrary feature pair
#' plotPairMulti(Vicari$TranscriptOutcomes, Vicari$MetaboliteOutcomes,
#'     Vicari$TranscriptCoords, Vicari$MetaboliteCoords,
#'     normX = "log", normY = "log", features = c("mt.Nd2", "X555.20713")
#' )
plotTopPair <- function(results, ..., normX = results$normX, normY = results$normY,
    topRank = 1, parameter = "Intercept") {
    stopifnot(is.numeric(topRank))
    if (!results$multi) {
        topFeats <- sund(rownames(results$result)[topRank])
        plotPairSingle(
            features = topFeats, assayX = results$assayX,
            assayY = results$assayY, normX = normX,
            normY = normY, ...
        )
    } else {
        stopifnot(parameter %in% names(results$result))
        topFeats <- sund(rownames(results$result[[parameter]])[topRank])
        plotPairMulti(
            features = topFeats, assayX = results$assayX,
            assayY = results$assayY, normX = normX,
            normY = normY, ...
        )
    }
}
#' @rdname plotTopPair
#' @export
#' @inheritParams plotPairSingle
#' @order 3
#' @param theme the ggplot2 theme
plotPairMulti <- function(
      Xl, Yl, Cxl, Eyl, features, normX = c("none", "rel", "log"),
      normY = c("none", "rel", "log"), size = 1.25, assayX, assayY, theme = theme_bw()
) {
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
    dfList <- Reduce(f = rbind, lapply(names(Xl), function(nam) {
        X <- normMat(Xl[[nam]], normX)
        Y <- normMat(Yl[[nam]], normY)
        coordMat <- rbind(Cxl[[nam]][rownames(X), ], Eyl[[nam]][rownames(Y), ])
        colnames(coordMat) <- c("x", "y")
        data.frame(
            "outcome" = c(
                scaleHelpFun(X, feat = features[1]),
                scaleHelpFun(Y, feat = features[2])
            ),
            "image" = nam, coordMat,
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
#' @param results Results returned by \link{sbivarSingle} or \link{extractResultsMulti}
#' @param x,y Outcome vectors
#' @param normX,normY Character strings, indicating what normalization is required
#' for X and Y matrices, respectively, before plotting, see details.
#' @param size Point size
#' @param features Feature vector of length 2 to be plotted
#' @rdname plotTopPair
#' @export
#' @order 2
plotPairSingle <- function(
      X, Y, Cx, Ey, features, normX = c("none", "rel", "log"),
      normY = c("none", "rel", "log"), assayX, assayY, ...
) {
    stopifnot(length(features) == 2)
    if (inherits(X, "SpatialExperiment")) {
        Cx <- SpatialExperiment::spatialCoords(X)
        X <- assayT(X, assayX)
    }
    if (inherits(Y, "SpatialExperiment")) {
        Ey <- SpatialExperiment::spatialCoords(Y)
        Y <- assayT(Y, assayY)
    }
    features <- make.names(features)
    colnames(X) <- make.names(colnames(X))
    colnames(Y) <- make.names(colnames(Y))
    foo <- checkInputSingle(X, Y, Cx, Ey)
    normX <- match.arg(normX)
    normY <- match.arg(normY)
    X <- normMat(X, normX)
    Y <- normMat(Y, normY)
    plotPairSingleVectors(
        x = scaleHelpFun(feat = features[1], X),
        y = scaleHelpFun(feat = features[2], Y),
        Cx = Cx[rownames(X), ], Ey = Ey[rownames(Y), ], modalityNames = features, ...
    )
}
#' @rdname plotTopPair
#' @param modalityNames Names to be given to the modalities,
#' appearing in the strip text of the columns. For plotTopPair() and
#' plotPairSingle(), the feature names are used.
#' @order 4
plotPairSingleVectors <- function(x, y, Cx, Ey, size = 1.25,
    modalityNames = c("Modality X", "Modality Y"), theme = theme_bw(), ...) {
    theme_set(theme)
    stopifnot(length(x) == nrow(Cx), length(y) == nrow(Ey), ncol(Ey) == 2, ncol(Cx) == 2)
    coordMat <- rbind(Cx, Ey)
    colnames(coordMat) <- c("x", "y")
    plotDf <- data.frame(
        "outcome" = c(x, y), coordMat,
        "feature" = rep(modalityNames, times = c(length(x), length(y)))
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
