#' Wrapper function to fit GAMs and find correlations and standard error for data lists
#'
#' @inheritParams wrapGAMs
#' @inheritParams sbivarMulti
#'
#' @returns A list named like Xl, containing all results
wrapGAMsMulti = function(Xl, Yl, Cxl, Eyl, families, n_points_grid){
    lapply(selfName(names(Xl)), function(nam){
        out <- wrapGAMs(Xl[[nam]], Yl[[nam]], Cxl[[nam]], Eyl[[nam]],
                 families = families, n_points_grid = n_points_grid)
        colnames(out)[seq_len(2)] = c("est", "se")
        return(out[,c("est", "se")])
    })
}
