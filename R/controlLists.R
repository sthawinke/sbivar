#' Create default argument lists
#'
#' Define default argument lists for the MoransISingle(), GAMsSingle() and GPsingle() functions
#'
#' @inheritParams MoransISingle
#' @order 1
#' @seealso \link{MoransISingle, fitGAM, fitGP}
#' @rdname controlLists
moranControl <- function(
      etas = c(5e-6, 2e-4, 2e-2), findMaxW = FALSE, returnSEsMoransI = TRUE,
      variogramModels = c("Exp", "Lin"), width = cutoff / 15, cutoff = sqrt(2) / 3,
      wo = c("Gauss", "nn"), numNNs = c(4, 8, 24)
) {
    list(
        etas = etas, findMaxW = findMaxW, returnSEsMoransI = returnSEsMoransI,
        variogramModels = match.arg(variogramModels), width = width, cutoff = cutoff,
        wo = match.arg(wo), numNNs = numNNs
    )
}
#' @inheritParams fitGAM
#' @order 2
#' @rdname controlLists
GAMcontrol <- function(Gamm = FALSE, Cx, n_points_grid = 6e2,
    families = list("X" = gaussian(), "Y" = gaussian())) {
    list(Gamm = Gamm, n_points_grid = n_points_grid, families = families)
}
#' @inheritParams fitGP
#' @order 3
#' @importFrom nlme lmeControl
#' @rdname controlLists
GPcontrol <- function(
      GPmethod = c("REML", "ML"), Cx, Quants = c(0.005, 0.5), numLscAlts = 5,
      optControl = lmeControl(
          opt = "optim", maxIter = 5e2, msMaxIter = 5e2,
          niterEM = 1e3, msMaxEval = 1e3
      )) {
    list(
        GPmethod = match.arg(GPmethod), Quants = Quants, numLscAlts = numLscAlts,
        optControl = optControl
    )
}
