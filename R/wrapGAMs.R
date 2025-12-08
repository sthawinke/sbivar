#' Wrapper function to fit GAMs and test for all possible combinations
#'
#' @inheritParams sbivarSingle
#' @param families A vector of length 2 giving the distributional families
#' for the outcome values. See details of \link{sbivarSingle}.
#' @param n_points_grid The number of points in the new grid for the GAMs to be
#' evaluated on.
#' @param multi Flag to indicate this is part of a multi-image analysis.
#' It only matters for the messages printed.
#' @returns A named list of results
wrapGAMs = function(X, Y, Cx, Ey, families, n_points_grid, verbose, multi = FALSE){
    if(verbose){
        message("Fitting GAMs for first modality (", ncol(X), " features)")
    }
    gamsx = fitManyGAMs(mat = X, coord = Cx, family = families[["X"]], modality = "X")
    if(verbose){
        message("Fitting GAMs for second modality (", ncol(Y), " features)")
    }
    gamsy = fitManyGAMs(mat = Y, coord = Ey, family = families[["Y"]], modality = "Y")
    ng = buildNewGrid(Cx = Cx, Ey = Ey, n_points_grid = n_points_grid)
    if(verbose){
        numTests = ncol(X)*ncol(Y)
        if(multi){
            message("Estimating all ", numTests,
                    " pairwise correlations on fitted GAMs ...")
        } else {
            message("Performing all ", numTests,
                    " pairwise tests on fitted GAMs ...")
            }
    }
    out = vapply(selfName(names(gamsx)), function(featx){
        predx <- vcovPredGam(gamsx[[featx]], newdata = ng)
        out = vapply(selfName(names(gamsy)), FUN.VALUE = double(3), function(featy){
            testGAM(predx = predx, modely = gamsy[[featy]], modelx = gamsx[[featx]],
                    predy = vcovPredGam(gamsy[[featy]], newdata = ng))
        })
        if(verbose)
            printProgress(featx, colnames(X))
        return(out)
    }, FUN.VALUE = matrix(0, nrow = 3, ncol = length(gamsy)))
    #Reformat to long format
    t(matrix(c(out), 3, length(gamsx)*length(gamsy),
             dimnames = list(c("corxy", "se.corxy", "pVal"),
              paste(rep(names(gamsx), each = length(gamsy)), rep(names(gamsy), times = length(gamsx)),
                    sep = "__"))))
}
