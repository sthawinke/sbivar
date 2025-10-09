#' Wrapper function to fit GAMs and find correlations and standard error for data lists
#'
#' @param families,n_points_grid See \link{wrapGAMs}
#' @inheritParams sbivarMulti
#'
#' @returns A list named like Xl, containing all results
wrapGAMsMulti = function(Xl, Yl, Cxl, Eyl, families, n_points_grid, verbose){
    lapply(selfName(names(Xl)), function(nam){
        if(verbose)
            printIteration(nam, names(Xl))
        out <- wrapGAMs(Xl[[nam]], Yl[[nam]], Cxl[[nam]], Eyl[[nam]],
                 families = families, n_points_grid = n_points_grid, verbose = verbose)
        colnames(out)[seq_len(2)] = c("est", "se")
        return(out[,c("est", "se")])
    })
}
