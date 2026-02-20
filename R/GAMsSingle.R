#' Fit univariate GAMs and test bivariate combinations
#'
#' Fit univariate GAMs with 2D cubic smoothing splines for both modalities, and test for correlation between all bivariate combinations.
#'
#' @inheritParams sbivarSingle
#' @param families A vector of length 2 giving the distributional families
#' for the outcome values. See details of \link{sbivarSingle}.
#' @param n_points_grid The number of points in the new grid for the GAMs to be
#' evaluated on.
#' @returns A named list of results
#' @inheritParams MoransISingle
GAMsSingle <- function(X, Y, Cx, Ey, families, n_points_grid, verbose, featuresX, featuresY, findVariances = TRUE) {
    if (verbose) {
        message("Fitting GAMs for first modality (", length(featuresX), " features) ...")
    }
    gamsx <- fitManyGAMs(mat = X, coord = Cx, family = families[["X"]], modality = "X", features = featuresX)
    if (verbose) {
        message("Fitting GAMs for second modality (", length(featuresY), " features) ...")
    }
    gamsy <- fitManyGAMs(mat = Y, coord = Ey, family = families[["Y"]], modality = "Y", features = featuresY)
    ng <- buildNewGrid(Cx = Cx, Ey = Ey, n_points_grid = n_points_grid)
    if (verbose) {
        numTests <- length(featuresX) * length(featuresY)
        message("Performing all ", numTests, " pairwise tests on fitted GAMs ...")
    }
    Nrow <- if (findVariances) 3 else 1
    out <- vapply(selfName(names(gamsx)), function(featx) {
        predx <- vcovPredGam(gamsx[[featx]], newdata = ng, findVariances = findVariances)
        out <- vapply(selfName(names(gamsy)), FUN.VALUE = double(Nrow), function(featy) {
            testGAM(
                predx = predx, modely = gamsy[[featy]], modelx = gamsx[[featx]],
                predy = vcovPredGam(gamsy[[featy]], newdata = ng, findVariances = findVariances),
                findVariances = findVariances
            )
        })
        printProgress(featx, featuresX, verbose)
        return(out)
    }, FUN.VALUE = matrix(0, nrow = Nrow, ncol = length(gamsy)))
    # Reformat to long format
    t(matrix(c(out), Nrow, length(gamsx) * length(gamsy),
        dimnames = list(
            c("corxy", if (findVariances) c("se.corxy", "pVal")),
            paste(rep(names(gamsx), each = length(gamsy)), rep(names(gamsy), times = length(gamsx)),
                sep = "__"
            )
        )
    ))
}
