#' @title Plot a feature pair consisting of a single molecule point pattern and a quantitative outcome
#' @description Plot a chosen feature pair, or the highest ranking feature pair,
#' for a single image
#' @inheritParams plotTopPair
#' @param modalityNames Names to be given to the modalities,
#' appearing in the strip text of the columns. For plotTopPairPPP() and
#' plotPairSinglePPP(), the feature names are used.
#' @param sideBySide Boolean indicating the patterns are to be plotted side by side, or on top of each other if FALSE
#' @param sizeX,sizeY Desired point sizes for corresponding modalities
#' @seealso \link{plotTopPair}, \link{sbivarSinglePPP}
#' @return A ggplot object
#' @order 3
#' @export
#' @rdname plotTopPairPPP
#' @note The point pattern is coloured blue, the value corresponding to 1
#' @examples
#' n <- 8e1
#' m <- 9e1
#' p <- 3
#' k <- 2
#' Y <- matrix(rnorm(m * k), m, k,
#'     dimnames =
#'         list(paste0("sampleY", seq_len(m)), paste0("Y", seq_len(k)))
#' )
#' Ey <- matrix(runif(m * 2), m, 2, dimnames = list(rownames(Y), c("x", "y")))
#' # Single image analysis on synthetic point pattern + quantitative modality
#' library(spatstat.random)
#' lambda <- 8e1
#' p <- 3
#' PPP <- rmpoispp(lambda, types = paste0("X", seq_len(p)), win = owin(c(0, 1), c(0, 1)))
#' marks(PPP, drop = FALSE) <- data.frame(feature = marks(PPP)) # Make sure marks is a dataframe
#' resMoransIppp <- sbivar(X = PPP, Y, Ey, method = "Moran's I")
#' plotTopPairPPP(resMoransIppp, X = PPP, Y = Y, Ey = Ey)
#' # For overlay, do:
#' plotPairPPP(X = PPP, Y = Y, Ey = Ey, features = c("X1", "Y1"), sideBySide = FALSE)
plotPairSinglePPvec <- function(Cx, y, Ey, sizeX = 0.005, sizeY = .01, sideBySide = TRUE,
    modalityNames = c("Modality X", "Modality Y"), theme = theme_bw(), ...) {
    theme_set(theme)
    stopifnot(length(y) == nrow(Ey), ncol(Ey) == 2, ncol(Cx) == 2)
    plotDfX <- data.frame(row.names = seq_len(nrow(Cx)),
        Cx, "outcome" = 1, "size" = sizeX,
        "feature" = modalityNames[1]
    )
    plotDfY <- data.frame(
        Ey, "outcome" = scaleZeroOne(y), "size" = sizeY,
        "feature" = modalityNames[2]
    )
    p <- if (sideBySide) {
        plotDf = rbind(plotDfX, plotDfY)
        plotDf$feature = factor(plotDf$feature, levels = modalityNames, ordered = TRUE)
        ggplot(data = plotDf, aes(x = x, y = y, col = outcome, size = size)) +
            geom_point() +
            facet_grid(~feature)
    } else {
        ggplot(data = plotDfY, aes(x = x, y = y, col = outcome, size = size)) +
            geom_point() +
            geom_point(inherit.aes = FALSE, data = plotDfX, aes(x = x, y = y, size = size))
    }
    p + scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("x coordinate") +
        ylab("y coordinate") +
        coord_fixed() +
        guides(size="none") +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
}
#' @rdname plotTopPairPPP
#' @order 2
#' @export
#' @importFrom spatstat.geom coords
plotPairPPP <- function(X, Y, Ey, features, normY = "none", ...) {
    Y <- normMat(Y, normY)
    plotPairSinglePPvec(Cx = coords(subset.ppp(X, marks(X, drop = FALSE)$feature == features[1])),
        y = Y[, features[2]],
        Ey = Ey, modalityNames = features, ...
    )
}
#' @export
#' @order 1
#' @rdname plotTopPairPPP
#' @importFrom spatstat.geom is.ppp
plotTopPairPPP <- function(result, X, Y, Ey, topRank = 1, ...) {
    stopifnot(is.numeric(topRank), topRank >= 1, is.ppp(X))
    plotPairPPP(X,
        Y = Y, Ey = Ey, normY = result$normY,
        features = unlist(result$result[topRank, c("Modality_X", "Modality_Y")]), ...
    )
}
