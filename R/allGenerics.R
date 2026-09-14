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
#' @examples
#' # Single image analysis on synthetic data
#' n <- 8e1
#' m <- 9e1
#' p <- 3
#' k <- 2
#' X <- matrix(rnorm(n * p), n, p,
#'     dimnames =
#'         list(paste0("sampleX", seq_len(n)), paste0("X", seq_len(p)))
#' )
#' Y <- matrix(rnorm(m * k), m, k,
#'     dimnames =
#'         list(paste0("sampleY", seq_len(m)), paste0("Y", seq_len(k)))
#' )
#' Cx <- matrix(runif(n * 2), n, 2, dimnames = list(rownames(X), c("x", "y")))
#' Ey <- matrix(runif(m * 2), m, 2, dimnames = list(rownames(Y), c("x", "y")))
#' resMoransI <- sbivar(X, Y, Cx, Ey, method = "Moran's I")
#' resGAMs <- sbivar(X, Y, Cx, Ey, method = "GAMs")
#' Y2 <- matrix(rnorm(n * k), n, k,
#'     dimnames =
#'         list(paste0("sampleX", seq_len(n)), paste0("Y", seq_len(k)))
#' )
#' resModtTestJoint <- sbivar(X, Y2, Cx, method = "Modified t-test")
#' # Single image analysis on synthetic data, converted to SpatialExperiment
#' library(SpatialExperiment)
#' seX <- SpatialExperiment(
#'     assays = list("transcripts" = t(X)),
#'     spatialCoords = Cx
#' )
#' seY <- SpatialExperiment(
#'     assays = list("metabolites" = t(Y)),
#'     spatialCoords = Ey
#' )
#' resModtGPs <- sbivar(seX, seY,
#'     assayX = "transcripts", assayY = "metabolites",
#'     method = "GPs"
#' )
setGeneric("sbivar", function(X, ...) standardGeneric("sbivar"))
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
#' @param parameter The linear model parameter used to find the feature with the strongest effect.
#' The default is the intercept, i.e. the overall effect.
#' @param topRank An integer, the feature pair with the rank-th smallest p-value is plotted
#' @param ... passed onto lower level functions
#' @param scaleBySampleSums A boolean, should the size of the spots be scaled by their sample sum, e.g. library size or total ion count? Recommended to reflect differences in certainty depending on sample sums.
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
#'     dimnames =
#'         list(paste0("sampleX", seq_len(n)), paste0("X", seq_len(p)))
#' )
#' Y <- matrix(rnorm(m * k), m, k,
#'     dimnames =
#'         list(paste0("sampleY", seq_len(m)), paste0("Y", seq_len(k)))
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
#'     Vicari$TranscriptCoords, Vicari$MetaboliteCoords,
#'     normX = "rel", normY = "rel", features = c("Gnas", "Tocopherol")
#' )
setGeneric(
    "plotTopPair",
    function(x, topRank = 1, parameter = "Intercept", scaleBySampleSums = FALSE,
             normX = x@normX, normY = x@normY, ...) standardGeneric("plotTopPair")
)
#' @export
#' @rdname plotGAMs
setGeneric(
    "plotGAMsTopPair",
    function(x, topRank = 1, parameter = "Intercept",
              ...) standardGeneric("plotGAMsTopPair")
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
#' @importFrom openxlsx createWorkbook writeData addWorksheet saveWorkbook getSheetNames
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
