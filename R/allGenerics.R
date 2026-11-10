#' Spatial bivariate association analysis
#' @description
#' Perform a bivariate spatial association analysis, either on a single or multiple images.
#' Depending on the input, the workhorse functions \link{sbivarSingle}
#' (single-image) or \link{sbivarMulti} (multi-image) are called.
#' @rdname sbivar
#' @importFrom methods setGeneric setMethod
#' @export
#' @param X the input object, containing measurements of the first modality, see methods('sbivar')
#' @param ... additional constructor and analysis arguments
#' @return A list containing the analysis result, along with parameters used in the analysis
#' @seealso \link{sbivarSingle}, \link{sbivarMulti}
#' @note Feature names are converted using the make.names() function
#' @examples
#' # Single image analysis on synthetic data
#' n <- 8e1
#' m <- 9e1
#' p <- 3
#' k <- 2
#' X <- matrix(rnorm(n * p), n, p,
#'   dimnames =
#'     list(paste0("sampleX", seq_len(n)), paste0("X", seq_len(p)))
#' )
#' Y <- matrix(rnorm(m * k), m, k,
#'   dimnames =
#'     list(paste0("sampleY", seq_len(m)), paste0("Y", seq_len(k)))
#' )
#' Cx <- matrix(runif(n * 2), n, 2, dimnames = list(rownames(X), c("x", "y")))
#' Ey <- matrix(runif(m * 2), m, 2, dimnames = list(rownames(Y), c("x", "y")))
#' resMoransI <- sbivar(X, Y, Cx, Ey, method = "Moran's I")
#' resGAMs <- sbivar(X, Y, Cx, Ey, method = "GAMs")
#' Y2 <- matrix(rnorm(n * k), n, k,
#'   dimnames =
#'     list(paste0("sampleX", seq_len(n)), paste0("Y", seq_len(k)))
#' )
#' resModtTestJoint <- sbivar(X, Y2, Cx, method = "Modified t-test")
#' # Single image analysis on synthetic data, converted to SpatialExperiment
#' library(SpatialExperiment)
#' seX <- SpatialExperiment(
#'   assays = list("transcripts" = t(X)),
#'   spatialCoords = Cx
#' )
#' seY <- SpatialExperiment(
#'   assays = list("metabolites" = t(Y)),
#'   spatialCoords = Ey
#' )
#' resModtGPs <- sbivar(seX, seY,
#'   assayX = "transcripts", assayY = "metabolites",
#'   method = "GPs"
#' )
setGeneric("sbivar", function(X, Y, ...) standardGeneric("sbivar"))
#' Return the top results found by sbivar
#' @param x The sbivarResults object
#'
#' @return A dataframe of top results
#' @export
setGeneric(
  "topPairs",
  function(x, ...) standardGeneric("topPairs")
)
#' @title Plot a the top feature pair according to a sbivar analysis
#' @description Plot a chosen feature pair, or the highest ranking feature pair,
#' for a single image or multiple images.
#' @inheritParams sbivarMulti
#' @inheritParams plotPairSingle
#' @inheritParams topPairs
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
#' The default is the intercept, i.e. the overall effect.
#' @param topRank An integer, the feature pair with the rank-th smallest p-value is plotted
#' @param ... passed onto lower level functions
#' @param scaleBySampleSums A boolean, should the size of the spots be scaled by their sample sum,
#' e.g. library size or total ion count? Recommended to reflect differences in certainty depending on sample sums.
#' @seealso \link{extractResultsMulti}, \link{fitLinModels}
#' @return A ggplot object
#' @export
#' @rdname plotTopPair
#' @details For sequence count data, such as transcriptomics, normalization
#' may be indicated to achieve clear plots (normX = "rel" or "log", see  \link{normMat}).
#' The normalization used for plotting is not necessarily the same as the one used for the analysis.
#' @examples
#' ### Single image
#' # Single image analysis on synthetic data
#' n <- 8e1
#' m <- 9e1
#' p <- 3
#' k <- 2
#' X <- matrix(rnorm(n * p), n, p,
#'   dimnames =
#'     list(paste0("sampleX", seq_len(n)), paste0("X", seq_len(p)))
#' )
#' Y <- matrix(rnorm(m * k), m, k,
#'   dimnames =
#'     list(paste0("sampleY", seq_len(m)), paste0("Y", seq_len(k)))
#' )
#' Cx <- matrix(runif(n * 2), n, 2, dimnames = list(rownames(X), c("x", "y")))
#' Ey <- matrix(runif(m * 2), m, 2, dimnames = list(rownames(Y), c("x", "y")))
#' resMoransI <- sbivar(X, Y, Cx, Ey, method = "Moran's I")
#' # Plot the feature pair with the most significant signal
#' plotTopPair(resMoransI, X, Y, Cx, Ey)
#' # Plot an arbitrary feature pair
#' plotPairSingle(X, Y, Cx, Ey, features = c("X1", "Y1"))
#' ### Multi image
#' data(Vicari)
#' # Plot an arbitrary feature pair
#' plotPairMulti(Vicari$TranscriptOutcomes, Vicari$MetaboliteOutcomes,
#'   Vicari$TranscriptCoords, Vicari$MetaboliteCoords,
#'   normX = "rel", normY = "rel", features = c("Gnas", "Tocopherol")
#' )
setGeneric(
  "plotTopPair",
  function(
    x, ..., topRank = 1, parameter = "Intercept", scaleBySampleSums = FALSE,
    normX = x@normX, normY = x@normY
  ) {
    standardGeneric("plotTopPair")
  }
)
#' Plot GAMs of top features
#' @description Fit GAMS for the top feature pair chosen, and plot the resulting spline surfaces
#' @export
#' @inheritParams plotTopPair
#' @inheritParams sbivarSingle
setGeneric(
  "plotGAMsTopPair",
  function(x, X, Y, Cx, Ey, topRank = 1, parameter = "Intercept", ...) {
    standardGeneric("plotGAMsTopPair")
  }
)
#' Write \emph{sbivar} results to an excel worksheet
#'
#' The results of single- or multi-image analysis are written to an excel spreadsheet
#' with separate tabs per parameter tested, sorted by increasing p-value.
#'
#' @param x The analysis results, from a call to \link{sbivar} (single-image)
#' or \link{extractResultsMulti} (multi-image)
#' @param file The file to write the results to
#' @param overwrite A boolean, should the file be overwritten if it exists already?
#' @param digits An integer, the number of significant digits to retain for the effect size,
#' raw and adjusted p-values
#' @param sigLevel The significance level threshold to use for the adjusted p-values,
#' only features exceeding the threshold are written to the file. Set this parameter to 1 to write all features
#' @seealso \link[openxlsx]{createWorkbook}
#' @details If no feature exceeds the significance threshold for a certain parameter,
#' an empty tab is created. For each fixed effect, a single tab is written.
#' The "baseline" tabs indicate the overall patterns, the other tabs are named after the fixed effects
#' and indicate departure from this baseline for this fixed effect.
#' @return Returns invisible with a message when writing operation successful,
#' otherwise throws a warning.
#' @export
#' @examples
#' example(sbivar, "sbivar")
#' # The significance level is set to 1 here for illustration,
#' # meaning that all feature pairs will be written to the spreadsheet.
#' # Single result
#' writeSbivarToXlsx(resGAMs, file = tmpFile <- tempfile(fileext = ".xlsx"), sigLevel = 1)
#' file.exists(tmpFile)
setGeneric(
  "writeSbivarToXlsx",
  function(x, file, overwrite = FALSE, digits = 3, sigLevel = 0.05) standardGeneric("writeSbivarToXlsx")
)
#' Fit linear models on result calculated for replicated images, and extract the results
#' @description Given result estimated from replicated images using \link{sbivarMulti},
#' fit linear models to determine significance.
#' To maintain interpretability of the intercept, continuous fixed effect variables are centered,
#' and sum coding is used for the categorical ones.
#'
#' @param result A result with the result of bivariate spatial association,
#' from a call to the \link{sbivar} function with multiple images
#' @param designDf A design dataframe
#' @param Formula A formula for the linear model to be fitted, can contain random effects.
#' @param Control A control list for lmerTest::lmer
#' @param inverseWeigh A boolean, should estimates be inverse weighed by variance for GAMs and Moran's I?
#' @param scaleByMax A boolean, should Moran's I be scaled by maximum values before plugging into the linear model?
#' @param verbose Should a message with number of linear models and cores be printed?
#' @returns For fitLinModels(), a list of linear models
#' @export
#' @details The left hand side of "Formula" can be provided or not,
#' but it will be overridden by "out" for downstream analysis
#'
#' @examples
#' # Multi-image analysis on Vicari data, using GAMs
#' data(Vicari)
#' # Subset to 5 images and 500 spots limit computing time
#' VicariMultiTest <- lapply(Vicari, function(x) lapply(x[1:5], function(y) y[1:500, ]))
#' VicariRes <- sbivar(VicariMultiTest$TranscriptOutcomes, VicariMultiTest$MetaboliteOutcomes,
#'   VicariMultiTest$TranscriptCoords, VicariMultiTest$MetaboliteCoords,
#'   method = "GAMs", control = list("families" = list("X" = mgcv::nb(), "Y" = Gamma(lin = "log")))
#' )
#' mouse <- substr(names(Vicari$TranscriptOutcomes)[1:5], 1, 10)
#' designDf <- data.frame("mouse" = mouse) # The design matrix
#' multiGAMLmms <- fitLinModels(VicariRes, designDf, Formula = ~ (1 | mouse))
#' # Extract the results
#' resGAMsMulti <- extractResultsMulti(multiGAMLmms, designDf = designDf)
#' topPairs(resGAMsMulti)
#' @importFrom lmerTest lmer
#' @importFrom stats formula terms model.matrix
#' @importFrom lme4 lmerControl .makeCC isSingular lFormula mkReTrms
#' @importFrom reformulas findbars nobars
#' @importFrom methods is
#' @importFrom smoppix centerNumeric named.contr.sum loadBalanceBplapply
#' @importFrom BiocParallel bplapply bpparam bpworkers
#' @seealso \link[lmerTest]{lmer}, \link[stats]{lm}, \link[sbivar]{sbivarMulti}, \link[stats]{p.adjust}
setGeneric(
  "fitLinModels",
  function(
    result, designDf, Formula, verbose = TRUE, inverseWeigh = FALSE, scaleByMax = TRUE,
    Control = lmerControl(
      check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL),
      check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
      check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)
    )
  ) {
    standardGeneric("fitLinModels")
  }
)
#' Extract results from linear models
#' @description extractResultsMulti() adds the results as matrix to the object, including adjusted p-values, and sorted by p-value.
#'
#' @param method Multiplicity correction method passed onto p.adjust
#'
#' @importFrom stats p.adjust anova
#' @importFrom methods is
#' @return The modified results object
#' @export
#' @rdname fitLinModels
#' @order 2
setGeneric(
    "extractResultsMulti",
    function(result, designDf, method = "BH"
    ) {
        standardGeneric("extractResultsMulti")
    }
)
