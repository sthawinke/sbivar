#' Fit GAMs and find correlations and standard error for data lists
#'
#' Wraps \link{GAMsSingle} for lists
#'
#' @param families,n_points_grid,includeGPsmooth,testSmooth See \link{GAMsSingle}
#' @inheritParams sbivarMulti
#' @inheritParams MoransIMulti
#'
#' @returns A list named like Xl, containing all results
GAMsMulti <- function(
      Xl, Yl, Cxl, Eyl, families, n_points_grid, verbose,
      includeGPsmooth, testSmooth, findVariances = FALSE
) {
    lapply(selfName(names(Xl)), function(nam) {
        if (verbose) {
            printIteration(nam, names(Xl))
        }
        out <- GAMsSingle(Xl[[nam]], Yl[[nam]], Cxl[[nam]], Eyl[[nam]],
            families = families, n_points_grid = n_points_grid, includeGPsmooth = includeGPsmooth,
            verbose = FALSE, findVariances = findVariances, featuresX = colnames(Xl[[nam]]),
            featuresY = colnames(Yl[[nam]]), testSmooth = testSmooth
        )
        return(list("res" = out[, c("corxy", if (findVariances) "se.corxy"), drop = FALSE]))
    })
}
