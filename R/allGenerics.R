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
#' resModtTestJoint <- sbivar(X, Y2, Cx, method = "Modified")
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
