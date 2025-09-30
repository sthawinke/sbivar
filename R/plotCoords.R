#' Plot the coordinates of two omics modalities
#'
#' Plot the coordinates of two modalities onto the same coordinate framework, in two different colours.
#' This is a useful check of the alignment and overlap.
#' @details plotCoordsMulti() is a wrapper for plotCoords for lists of coordinates,
#' and requires the user to set par(mar = ) appropriately, so all plots are shown.
#'
#' @inheritParams sbivarSingle
#' @param cex Expansion factor
#' @param pch,pchY Point shapes for x and y
#' @param ... passed onto plot()
#'
#' @returns Plots to the plotting device
#' @export
#'
#' @examples
#' example(sbivarSingle, "sbivar")
#' plotCoords(Cx, Ey)
#' #For multiple coordinates
#' example(sbivarMulti, "sbivar")
#' par(mfrow = c(2,3))
#' foo = lapply(names(Cxl), function(nam){
#'     plotCoords(Cxl[[nam]], Eyl[[nam]], main = nam)
#' })
#' par(mfrow=c(1,1))
#' @importFrom graphics points
plotCoords = function(Cx, Ey, pch = 1, pchY = 3, cex = 0.8, ...){
    plot(x = as.matrix(Cx), xlim = range(c(Cx[,1], Ey[,1])), ylim = range(c(Cx[,2], Ey[,2])),
         asp = 1, pch = pch, cex = cex, xlab = "x", ylab = "y", ...)
    points(as.matrix(Ey), col = "blue", pch = pchY, cex = cex)
}
#' @inheritParams sbivarMulti
#' @export
#' @rdname plotCoords
plotCoordsMulti = function(Cxl, Eyl, ...){
    stopifnot(identical(names(Cxl),names(Eyl)))
    foo = lapply(names(Cxl), function(nam) {
        plotCoords(Cxl[[nam]], Eyl[[nam]], main =nam, ...)
    })
    return(invisible())
}
