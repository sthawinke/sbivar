#' Wrapper function to fit GAMs and test for all possible combinations
#'
#' @inheritParams sbivarSingle
#'
#' @returns A named list of results
#' @importFrom smoppix loadBalanceBplapply
#' @importFrom BiocParallel bplapply
wrapGAMs = function(X, Y, Cx, Ey, families, n_points_grid){
    gamsx = fitManyGAMs(mat = X, coord = Cx, family = families[["X"]])
    gamsy = fitManyGAMs(mat = Y, coord = Ey, family = families[["Y"]])
    ng = buildNewGrid(Cx = Cx, Ey = Ey, n_points_grid = n_points_grid)
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
                             makeNames(names(gamsx), names(gamsy)))))
}
