#' Calculate ibivariate Moran's I between two modality matrix, with variance and p-value
#'
#' The variance calculation requires estimation of the spatial autocorrelation structure of every feature separately
#'
#' @inheritParams sbivarSingle
#'
#' @returns A dataframe of results sorted by p-value, also containing the estimated Moran's I statistic and its variance.
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
wrapMoransI = function(X, Y, Cx, Ey, wo, verbose, ...){
    n = nrow(X);m = nrow(Y)
    featGrid = expand.grid("featX" = colnames(X), "featY" = colnames(Y))
    if(verbose){
        message("Performing tests on bivariate Moran's I for ", ncol(X)*ncol(Y), " feature pairs")
    }
    variogramsX = matheronVariograms(X, Cx, ...)
    variogramsY = matheronVariograms(Y, Ey, ...)
    out = simplify2array(loadBalanceBplapply(seq_len(nrow(featGrid)), function(i){
        unlist(modified.ttest(X[, featGrid[i, "featX"]], Y[, featGrid[i, "featY"]],
                              Cx)[c("corr", "ESS", "p.value")])
    }))
    colnames(out) = makeNames(colnames(X), colnames(Y))
    rownames(out) = c("Correlation", "Effective sample size", "pVal")
    t(out)
}
#' @importFrom gstat gstat variogram
matheronVariograms <- function(X, Cx, cutoff, width) {
    # coords: matrix or data.frame of coordinates (n x 2)
    # values: numeric vector of observations (length n)
    # cutoff: maximum distance to consider (optional)
    # width: bin width for distance classes (optional)
    df = data.frame(x = coords[,1], y = coords[,2], z = X)
    sp::coordinates(df) <- ~x + y
    g <- NULL
    for (nm in colnames(X)) {
        g <- gstat(g, id = nm, formula = as.formula(paste0(nm,"~1")), data=df)
    }
    # Compute empirical semivariogram using Matheron’s estimator
    v <- variogram(g, cutoff = cutoff, width = width, cross = FALSE)
    return(v)
}
