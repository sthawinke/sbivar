#' Estimate measures of bivariate spatial association for multiple images
#'
#' This function calculates measures of spatial association for every image.
#' The resulting estimates can then be analysed further using the \link{fitLinModels} function.
#'
#' @param Xl,Yl Lists of matrices of omics measurements
#' @param Cxl,Eyl Lists of corresponding coordinate matrices of dimension two
#' @param findMaxW A boolean, should maximum values be found for bivariate Moran's I
#' @param findVariances A boolean, should variances be found for bivariate Moran's I and GAMs
#' @returns A \link{sbivarResults-class} object with multi = TRUE
#' @note All methods use multithreading on the cluster provided using the BiocParallel package
#' @inheritParams sbivarSingle
#' @importFrom BiocParallel bpparam bpworkers
#' @seealso \link{fitLinModels}, \link{MoransIMulti}, \link{correlationsMulti}, \link{GAMsMulti}
sbivarMulti <- function(Xl, Yl, Cxl, Eyl, method = c("Moran's I", "GAMs", "Correlation"),
                        control, featuresX = getFeaturesList(Xl), featuresY = getFeaturesList(Yl),
                        normX = c("none", "rel", "log"), normY = c("none", "rel", "log"), pseudoCount = 1e-8,
                        verbose = TRUE, findVariances = FALSE, findMaxW = TRUE) {
  method <- match.arg(method)
  normX <- match.arg(normX)
  normY <- match.arg(normY)
  controlvals <- switch(method,
    "Moran's I" = moranControl(),
    "GAMs" = GAMcontrol()
  )
  if (!missing(control)) {
    controlvals[names(control)] <- control
  }
  invisible(checkControlVals(method, controlvals))
  stopifnot(
    is.logical(verbose), is.character(method), is.character(featuresX), is.character(featuresY)
  )
  ccs <- (method != "Correlation")
  invisible(checkInputMulti(Xl, Yl, Cxl, Eyl, checkCoords = ccs))
  if (ccs) {
    if (missing(Eyl)) {
      message("Only one coordinate matrix list Cxl supplied, and dimensions of X and Y do match.
                 Performing an analysis with joint coordinate sets.")
      Eyl <- Cxl
    }
    Cxl <- lapply(Cxl, tmpFun <- function(cx) {
      colnames(cx) <- c("x", "y")
      cx
    })
    Eyl <- lapply(Eyl, tmpFun)
  } else if (!missing(Cxl)) {
    warning("Correlation analysis will ignore coordinate matrices provided.
                Consider providing another 'method' argument for a full spatial analysis", immediate. = TRUE)
  }
  # Normalization
  Xl <- lapply(Xl, normMat, normX, pseudoCount, warnDrop = FALSE)
  Yl <- lapply(Yl, normMat, normY, pseudoCount, warnDrop = FALSE)
  if (ccs) {
    Cxl <- lapply(selfName(names(Cxl)), function(i) {
      Cxl[[i]][rownames(Xl[[i]]), ]
    })
    Eyl <- lapply(selfName(names(Eyl)), function(i) {
      Eyl[[i]][rownames(Yl[[i]]), ]
    })
  }
  if (verbose) {
    message(
      "Starting sbivar analysis (", method, ") of ", length(Xl), " images on ",
      bpworkers(bpparam()), " computing cores"
    )
  }
  res <- if (method == "Moran's I") {
    wo <- match.arg(controlvals$wo, choices = formals(moranControl)$wo)
    controlvals$findMaxW <- findMaxW
    if (wo == "nn" & !requireNamespace("RANN", quietly = TRUE)) {
      stop("Install 'RANN' package for nearest neighbour weighting matrices")
    }
    MoransIMulti(Xl, Yl, Cxl, Eyl,
      featuresX = featuresX, featuresY = featuresY,
      control = controlvals, verbose = verbose, findVariances = findVariances
    )
  } else if (method == "GAMs") {
    GAMsMulti(Xl, Yl, Cxl, Eyl,
      featuresX = featuresX, featuresY = featuresY,
      findVariances = findVariances,
      control = controlvals, verbose = verbose
    )
  } else if (method == "Correlation") {
    correlationsMulti(Xl, Yl,
      featuresX = featuresX, featuresY = featuresY, verbose = verbose
    )
  }
  out <- if (method == "Moran's I") {
    wParams <- switch(controlvals$wo,
      "Gauss" = controlvals$etas,
      "nn" = controlvals$numNNs
    )
    new("sbivarResultsMoransI",
      "result" = lapply(res, function(x) x$res), "method" = method, "normX" = normX, "multi" = TRUE,
      "normY" = normY, "maxIxy" = if (controlvals$findMaxW) vapply(res, FUN.VALUE = wParams, function(x) x$maxIxy),
      "wo" = controlvals$wo, "returnSEsMoransI" = findVariances,
      "wParams" = wParams
    )
  } else if (method == "GAMs") {
    new("sbivarResultsGAMs",
      "result" = res, "method" = method, "multi" = TRUE,
      "normX" = normX, "normY" = normY, "families" = controlvals$families
    )
  } else if (method == "Correlation") {
    new("sbivarResults",
      "result" = res, "method" = method, "multi" = TRUE,
      "normX" = normX, "normY" = normY
    )
  }
  return(out)
}
