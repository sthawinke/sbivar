#' Wrapper function to fit GAMs and test for all possible combinations
#'
#' @inheritParams sbivarSingle
#'
#' @returns see \link{testManyGAMs}
wrapGAMs = function(X, Y, Cx, Ey, families, n_points_grid){
    xGams = fitManyGAMs(mat = X, coord = Cx, family = families[["X"]])
    yGams = fitManyGAMs(mat = Y, coord = Ey, family = families[["Y"]])
    ng = buildNewGrid(Cx = Cx, Ey = Ey, n_points_grid = n_points_grid)
    testManyGAMs(xGams, yGams, newGrid = ng)
}
