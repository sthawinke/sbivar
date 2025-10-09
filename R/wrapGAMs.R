#' Wrapper function to fit GAMs and test for all possible combinations
#'
#' @inheritParams sbivarSingle
#' @param families A vector of length 2 giving the distributional families
#' for the outcome values. See details of \link{sbivarSingle}.
#' @param n_points_grid The number of points in the new grid for the GAMs to be
#' evaluated on.
#' @returns A named list of results
#' @importFrom smoppix loadBalanceBplapply
#' @importFrom BiocParallel bplapply
wrapGAMs = function(X, Y, Cx, Ey, families, n_points_grid, verbose){
    if(verbose){
        message("Fitting GAMs for first modality (", ncol(X), " features)")
    }
    gamsx = fitManyGAMs(mat = X, coord = Cx, family = families[["X"]], modality = "X")
    if(verbose){
        message("Fitting GAMs for first modality (", ncol(Y), " features)")
    }
    gamsy = fitManyGAMs(mat = Y, coord = Ey, family = families[["Y"]], modality = "Y")
    ng = buildNewGrid(Cx = Cx, Ey = Ey, n_points_grid = n_points_grid)
    if(verbose){
        message("Performing all ", ncol(X)*ncol(Y),
                " pairwise tests on fitted GAMs ...")
    }
    out = loadBalanceBplapply(selfName(names(gamsx)), function(featx){
        predx <- vcovPredGam(gamsx[[featx]], newdata = ng)
        vapply(selfName(names(gamsy)), FUN.VALUE = double(3), function(featy){
            testGAM(predx = predx, modely = gamsy[[featy]], modelx = gamsx[[featx]],
                    predy = vcovPredGam(gamsy[[featy]], newdata = ng))
        })
    })
    #Reformat to long format
    t(matrix(unlist(out), 3, length(gamsx)*length(gamsy),
             dimnames = list(c("corxy", "se.corxy", "pVal"),
              paste(rep(names(gamsx), each = length(gamsy)), rep(names(gamsy), times = length(gamsx)),
                    sep = "__"))))
}
