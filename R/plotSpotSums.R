#' Plot the spot-wise sums of the two modalities across the slices
#'
#' This is an important exploratory plot to visualize spot-wise sums
#' (library sizes, total ion counts etc.) spatially
#'
#' @inheritParams sbivarMulti
#' @param viewNames The names of the views, showed in the column names
#' @param pointSize The point size for geom_point()
#' @param stripTextSize The size of the text in the strips (grey panels)
#'
#' @returns A ggplot2 object
#' @export
#'
#' @examples
#' data(Vicari)
#' plotSpotSums(
#'     Vicari$TranscriptOutcomes, Vicari$MetaboliteOutcomes,
#'     Vicari$TranscriptCoords, Vicari$MetaboliteCoords
#' )
plotSpotSums <- function(Xl, Yl, Cxl, Eyl, viewNames = c("Transcriptomics", "Metabolomics"),
    pointSize = 0.01, stripTextSize = 8.5) {
    baa <- checkInputMulti(Xl, Yl, Cxl, Eyl)
    foo <- base::do.call(rbind, lapply(names(Xl), function(x) {
        dfStx <- data.frame(Cxl[[x]], "log10_Library_size" = log10(rowSums(Xl[[x]])), what = "Transcriptomics")
        dfMsi <- data.frame(Eyl[[x]], "log10_Library_size" = log10(rowSums(Yl[[x]])), what = "Metabolomics")
        df <- rbind(dfStx, dfMsi)
        df$sample <- x
        df
    }))
    foo$what <- factor(foo$what, levels = viewNames)
    p <- ggplot(data = foo, aes(x = x, y = y, col = log10_Library_size)) +
        geom_point(size = pointSize) +
        facet_grid(sample ~ what, space = "free", scales = "free") +
        theme(strip.text = element_text(size = stripTextSize)) +
        xlab("Coordinate 1") +
        ylab("Coordinate 2") +
        scale_colour_gradient(
            low = "yellow", high = "blue",
            name = "Log10(sample sum)"
        ) +
        coord_fixed()
    return(p)
}
