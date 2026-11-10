#' Test for bivariate spatial association in a single image
#'
#' The tests can be applied to either joint or disjoint coordinate sets.
#'
#' @details If only Cx is supplied and X and Y have the same number of rows, a joint analysis is performed
#' If Cx and Ey are provided, and X and Y have the same number of rows, equality of Cx and Ey is checked.
#' If true, a joint analysis is run, with a warning.
#'
#' @param X,Y Matrices of omics measurements
#' @param Cx,Ey Coordinate matrices of dimension two, belonging to X and Y respectively
#' @param method A character string, indicating which method to apply
#' @param correlation Correlation structure, passed onto \link{fitGP} or \link{fitGAM}
#' @param verbose Should info on type of analysis be printed?
#' @param control Control arguments passed onto \link{MoransIsingle}, \link{fitGPs} or \link{fitGAMs}, replacing defaults by \link{MoransIsingle}, \link{fitGPs} or \link{fitGAMs}, respectively
#' @param normX,normY,pseudoCount Normalization parameters, passed onto \link{normMat}
#' @param featuresX,featuresY Features to be tested. Defaults to all features,
#' but specifying them allows to test a limited feature set,
#' while using the whole matrix to calculate library sizes as offset or for normalization.
#'
#' @details X and Y need to have rownames for matching to the coordinates,
#' and column names for identifying the features.
#' Cx and Ey must have rownames matching those in X and Y, and have two columns.
#' For GAMs, usually no normalization is needed, as the non-gaussianity is taken care of by
#' the outcome distribution, offset and link functions. Currently, identity,
#' inverse and log-link are implemented.
#'
#' @returns A \link{sbivarResults} object with multi = FALSE
#' @importFrom stats p.adjust
#' @importFrom methods is new
#' @importFrom nlme corRatio corGaus corSpher corExp corLin lmeControl
#' @importFrom BiocParallel bpparam bpworkers
#' @note All methods use multithreading on the cluster provided using the BiocParallel package.
#' @seealso \link{MoransISingle}, \link{ModTtestSingle}, \link{GAMsSingle}, \link{GPsSingle}
sbivarSingle <- function(X, Y, Cx, Ey, method = c("Moran's I", "GAMs", "Modified t-test", "GPs"),
                         normX = c("none", "rel", "log"), normY = c("none", "rel", "log"), pseudoCount = 1e-8,
                         correlation = corGaus(form = ~ x + y, nugget = TRUE, value = c(0.9 * max(apply(Cx, 2, function(x) diff(range(x)))), 0.25)), control, gpParams,
                         featuresX = colnames(X), featuresY = colnames(Y),
                         verbose = TRUE) {
  method <- match.arg(method)
  normX <- match.arg(normX)
  normY <- match.arg(normY)
  controlvals <- switch(method,
    "Moran's I" = moranControl(),
    "GPs" = GPcontrol(),
    "GAMs" = GAMcontrol()
  )
  if (!missing(control)) {
    controlvals[names(control)] <- control
  }
  variogramModels <- match.arg(controlvals$variogramModels,
    choices = formals(moranControl)$variogramModels,
    several.ok = TRUE
  )
  stopifnot(
    ncol(Cx) == 2, !is.null(colnames(Y)), all(featuresX %in% colnames(X)),
    all(featuresY %in% colnames(Y)), !anyDuplicated(featuresX), !anyDuplicated(featuresY),
    !is.null(colnames(X)), is.logical(verbose)
  )
  invisible(checkControlVals(method, controlvals, gpParams))
  if (method == "Moran's I") {
    wo <- match.arg(controlvals$wo, choices = formals(moranControl)$wo)
  } else if (method == "GPs") {
    GPmethod <- match.arg(controlvals$GPmethod, choices = formals(GPcontrol)$GPmethod)
  }
  invisible(checkInputSingle(X, Y, Cx, Ey))
  colnames(Cx) <- c("x", "y")
  X <- normMat(X, normX, pseudoCount)
  Y <- normMat(Y, normY, pseudoCount)
  Cx <- Cx[rownames(X), ] # Match coordinates to observations through rownames
  featuresX <- make.names(featuresX)
  featuresY <- make.names(featuresY)
  if (missing(Ey)) {
    if (nrow(X) != nrow(Y)) {
      stop(strwrapPaste("Only one coordinate matrix supplied, but sample size of X and Y do not match!
                 Please supply the 'Ey' argument and consider 'Moran's I' as method."))
    }
    # Run a joint analysis
    if (verbose) {
      message(strwrapPaste("Only one coordinate matrix Cx supplied, and dimensions of X and Y do match.
                 Performing an analysis with joint coordinate sets."))
    }
    Ey <- Cx
  } else if (method == "Modified t-test") {
    stop(strwrapPaste("Two coordinate matrices supplied, whereas the modified t-test requires only one.
            If coordinates are shared, omit Ey. If coordinates are disjoint, consider using method = 'MoransI'."))
  } else {
    Ey <- Ey[rownames(Y), ]
    colnames(Ey) <- c("x", "y")
  }
  if (((normX %in% c("rel", "log")) || (normY %in% c("rel", "log"))) && method == "GAMs") {
    warning(strwrapPaste("Normalizing data is not recommended for GAMs!
                try accounting for non-normality through the 'families' argument."), immediate. = TRUE)
  }
  if (verbose) {
    message(
      strwrapPaste(
        "Starting sbivar analysis of a single image on ",
        bpworkers(bpparam()), " computing cores"
      )
    )
  }
  out <- if (method == "Moran's I") {
    (moranRes <- MoransISingle(
      X = X, Y = Y, Cx = Cx, Ey = Ey, verbose = verbose, wo = controlvals$wo, numNNs = selfName(controlvals$numNNs),
      variogramModels = controlvals$variogramModels, etas = selfName(controlvals$etas), width = controlvals$width,
      returnSEsMoransI = controlvals$returnSEsMoransI,
      cutoff = controlvals$cutoff, findMaxW = controlvals$findMaxW,
      featuresX = featuresX, featuresY = featuresY
    ))$res
  } else if (method == "GAMs") {
    GAMsSingle(
      X = X, Y = Y, Cx = Cx, Ey = Ey, families = controlvals$families, n_points_grid = controlvals$n_points_grid,
      verbose = verbose, featuresX = featuresX, featuresY = featuresY, Gamm = controlvals$Gamm, correlation = correlation
    )
  } else if (method == "GPs") {
    GPsSingle(
      X = X, Y = Y, Cx = Cx, Ey = Ey, gpParams = gpParams, Quants = controlvals$Quants,
      GPmethod = controlvals$GPmethod, correlation = correlation, optControl = controlvals$optControl,
      numLscAlts = controlvals$numLscAlts, verbose = verbose, featuresX = featuresX, featuresY = featuresY
    )
  } else if (method == "Modified t-test") {
    sharedNames <- intersect(rownames(X), rownames(Y))
    ModTtestSingle(
      X = X[sharedNames, featuresX, drop = FALSE],
      Y = Y[sharedNames, featuresY, drop = FALSE],
      Cx = Cx[sharedNames, ], verbose = verbose
    )
  }
  out <- cbind(out, "pAdj" = p.adjust(out[, "pVal"], method = "BH"))
  out <- addFeatureColumn(out[order(out[, "pVal"]), , drop = FALSE])
  res <- if (method == "Moran's I") {
    new("sbivarResultsMoransI",
      "result" = out, "method" = method, "normX" = normX, "multi" = FALSE,
      "normY" = normY, "maxIxy" = moranRes$maxIxy, "wo" = controlvals$wo,
      "returnSEsMoransI" = controlvals$returnSEsMoransI,
      "wParams" = switch(controlvals$wo,
        "Gauss" = controlvals$etas,
        "nn" = controlvals$numNNs
      )
    )
  } else if (method == "GAMs") {
    new("sbivarResultsGAMs",
      "result" = out, "method" = method, "multi" = FALSE,
      "normX" = normX, "normY" = normY, "families" = controlvals$families,
      "correlation" = if (controlvals$Gamm) controlvals$correlation, "Gamm" = controlvals$Gamm
    )
  } else {
    new("sbivarResults",
      "result" = out, "method" = method, "multi" = FALSE,
      "normX" = normX, "normY" = normY
    )
  }
  return(res)
}
