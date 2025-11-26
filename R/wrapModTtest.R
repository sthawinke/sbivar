#' Perform modified t-tests for all pairs
#'
#' Then modified t-test is applied to all pairs, which tests for the significance
#' of the Pearson correlation while accounting for spatial autocorrelation \insertCite{Clifford1989}{sbivar}.
#'
#' @inheritParams sbivarSingle
#'
#' @returns A dataframe of results sorted by p-value, also containing effective sample size (ESS) and correlation estimate.
#' @importFrom SpatialPack modified.ttest
#' @seealso \link[SpatialPack]{modified.ttest}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
wrapModTtest = function(X, Y, Cx, verbose){
    n = nrow(X);m = nrow(Y)
   featGrid = expand.grid("featX" = colnames(X), "featY" = colnames(Y))
    if(verbose){
        message("Performing ", nrow(featGrid), " modified t-tests")
    }
    out = simplify2array(loadBalanceBplapply(seq_len(nrow(featGrid)), function(i){
        unlist(modified.ttest(X[, featGrid[i, "featX"]], Y[, featGrid[i, "featY"]],
                       Cx)[c("corr", "ESS", "p.value")])
    }))
    colnames(out) = makeNames(colnames(X), colnames(Y))
    rownames(out) = c("Correlation", "Effective sample size", "pVal")
    t(out)
}
