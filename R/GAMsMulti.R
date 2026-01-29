#' Fit GAMs and find correlations and standard error for data lists
#'
#' Wraps \link{GAMsSingle} for lists
#'
#' @param families,n_points_grid See \link{GAMsSingle}
#' @inheritParams sbivarMulti
#' @inheritParams MoransIMulti
#'
#' @returns A list named like Xl, containing all results
GAMsMulti <- function(Xl, Yl, Cxl, Eyl, families, n_points_grid, verbose, findVariances = FALSE) {
    lapply(selfName(names(Xl)), function(nam) {
        if (verbose) {
            printIteration(nam, names(Xl))
        }
        out <- GAMsSingle(Xl[[nam]], Yl[[nam]], Cxl[[nam]], Eyl[[nam]],
            families = families, n_points_grid = n_points_grid,
            verbose = FALSE, findVariances = findVariances
        )
        return(list("res" = out[, c("corxy", if (findVariances) "se.corxy"), drop = FALSE]))
    })
}
