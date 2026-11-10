#' Find all Moran's I statistics for a list of matrices
#'
#' Find all Moran's I values for lists of observations matrices from different modalities,
#' by calling the \link{MoransISingle} function.
#'
#' @inheritParams sbivarMulti
#' @inheritParams buildWeightMat
#' @inheritParams MoransISingle
#' @param ... passed onto \link{MoransISingle}
#'
#' @returns A list of Moran's I estimates, standard errors and maximum values
#' @seealso \link{MoransISingle}
MoransIMulti <- function(Xl, Yl, Cxl, Eyl, verbose, featuresX, featuresY, findVariances = FALSE, control, ...) {
  lapply(selfName(names(Xl)), function(nam) {
    if (verbose) {
      printIteration(nam, names(Xl))
    }
    MoransISingle(
      X = Xl[[nam]], Y = Yl[[nam]], Cx = Cxl[[nam]], Ey = Eyl[[nam]],
      verbose = FALSE, wo = control$wo, numNNs = selfName(control$numNNs),
      variogramModels = control$variogramModels, etas = selfName(control$etas), width = control$width,
      cutoff = control$cutoff, findMaxW = control$findMaxW, findVariances = findVariances,
      returnSEsMoransI = findVariances, featuresX = intersect(featuresX, colnames(Xl[[nam]])),
      featuresY = intersect(featuresY, colnames(Yl[[nam]])), ...
    )[c("res", "maxIxy")]
  })
}
