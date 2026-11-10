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
MoransIMulti <- function(Xl, Yl, Cxl, Eyl, verbose, featuresX, featuresY, control, ...) {
  lapply(selfName(names(Xl)), function(nam) {
    if (verbose) {
      printIteration(nam, names(Xl))
    }
    MoransISingle(
      X = Xl[[nam]], Y = Yl[[nam]], Cx = Cxl[[nam]], Ey = Eyl[[nam]],
      verbose = FALSE, wo = controlvals$wo, numNNs = selfName(controlvals$numNNs),
      variogramModels = controlvals$variogramModels, etas = selfName(controlvals$etas), width = controlvals$width,
      cutoff = controlvals$cutoff, findMaxW = controlvals$findMaxW, findVariances = findVariances,
      returnSEsMoransI = findVariances, featuresX = intersect(featuresX, colnames(Xl[[nam]])),
      featuresY = intersect(featuresY, colnames(Yl[[nam]])), ...
    )[c("res", "maxIxy")]
  })
}
