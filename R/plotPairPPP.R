#' @title Plot a feature pair consisting of a single molecule point pattern and a quantitative outcome
#' @description Plot a chosen feature pair, or the highest ranking feature pair,
#' for a single image
#' @inheritParams plotTopPair
#' @seealso \link{plotTopPair}, \link{sbivarSinglePPP}
#' @return A ggplot object
#' @order 1
#' @export
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
#' PPP <- rmpoispp(lambda, types = paste0("gene", seq_len(p)), win = owin(c(0, 1), c(0, 1)))
#' marks(PPP, drop = FALSE) <- data.frame(feature = marks(PPP)) # Make sure marks is a dataframe
#' resMoransIppp <- sbivar(X = PPP, Y, Ey, method = "Moran's I")
#' plotTopPairSinglePPP(resMoransIppp, PPP = X, Y = Y, Ey = Ey)
plotPairSinglePPvec <- function(Cx, y, Ey, sizeX = 0.1, sizeY = 1, shapeX = 15, shapeY = 1,
    modalityNames = c("Modality X", "Modality Y"), theme = theme_bw(), ...) {
    theme_set(theme)
    stopifnot(length(y) == nrow(Ey), ncol(Ey) == 2, ncol(Cx) == 2)
    plotDfY <- data.frame(
        "outcome" = scaleZeroOne(y), Ey,
        "feature" = modalityNames[2]
    )
    plotDfX <- data.frame(
        Cx,
        "feature" = modalityNames[1]
    )
    ggplot(data = plotDfY, aes(x = x, y = y, col = outcome)) +
        geom_point(size = size, shape = shapeY) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("x coordinate") +
        ylab("y coordinate") +
        coord_fixed() +
        geom_point(inherit.aes = FALSE, data = plotDfX, aes(x = x, y = y), size = sizeX, shape = shapeX)
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}
#' @importFrom spatstat.geom coords
plotPairPPP <- function(PPP, ...) {
    plotPairSinglePPmat(coords(PPP), ...)
}
#'@export
plotTopPairPPP <- function(result, X, Y, Ey, topRank = 1, ...) {
    stopifnot(is.numeric(topRank), topRank >= 1)
    Y <- normMat(Y, result$normY)
    plotPairSinglePPvec(coords(subset.ppp(X, "feature" == result$result[topRank, "Modality_X"])),
        y = Y[, result$result[topRank, "Modality_Y"]], Ey = Ey, ...
    )
}
