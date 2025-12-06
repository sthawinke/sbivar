#' Visualize different weighting functions
#' @description
#' The weighting functions resulting from different choices of decay parameters eta are visualized in a lineplot
#'
#' @param etas A vector of positive decay parameters
#' @param dists A set of distances, smaller than the square root of 2 since all coordinates are scaled to the unit square
#' @param palette The colour palette
#'
#' @returns Plots the weighting functions
#' @export
#'
#' @examples
#' exploreWeights(10^c(-5, -4, -3))
#' @importFrom grDevices palette
#' @importFrom graphics legend lines
exploreWeights = function(etas, dists = seq(0, 0.1, length.out = 1e3), palette = "paired"){
    stopifnot(all(dists<=sqrt(2)), all(etas>0))
    Pal = palette("paired")
    if(length(etas)>length(Pal)){
        stop("Insufficient colours in palette ", palette, "for supplied etas!")
    }
    plot(type = "n", x = 0, y = 0, xlim = range(dists), ylim = c(0,1), xlab = "distance", ylab = "weight")
    foo = lapply(seq_along(etas), function(i){
        lines(dists, exp(-dists^2/etas[i]), col = Pal[i])
    })
    legend("topright", legend = signif(etas, 3), lty =1, col = Pal[seq_along(etas)],
            title = "Decay parameter")
}
