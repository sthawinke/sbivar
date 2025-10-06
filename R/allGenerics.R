#' Perform a spatial bivariate association analysis
#' @description
#' Perform the spatial analysis of choice, either on a single or multiple images.
#' Depending on the input, the workhorse functions \link{sbivarSingle}
#' (single-image) or \link{sbivarMulti} (multi-image) are called.
#' @rdname sbivar
#' @importFrom methods setGeneric setMethod
#' @export
#' @param X the input object, see methods('sbivar')
#' @param ... additional constructor and analysis arguments
#' @return The analysis result
#' @seealso \link{sbivarSingle}, \link{sbivarMulti}
#' @examples
#' # Single image analysis on synthetic data
#' n=7e1;m=9e1;p=3;k=4
#' X = matrix(rnorm(n*p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
#' Y = matrix(rnorm(m*k), m, k, dimnames = list(NULL, paste0("Y", seq_len(k))))
#' Cx = matrix(runif(n*2), n, 2)
#' Ey = matrix(runif(m*2), m, 2)
#' colnames(Cx) = colnames(Ey) = c("x", "y")
#' resGAMs = sbivar(X, Y, Cx, Ey, method = "GAMs")
#' resModtTest = sbivar(X, Y, Cx, Ey, method = "Modified")
#' resModtTestJoint = sbivar(X, Y[seq_len(nrow(X)),], Cx, method = "Modified")
#' # Single image analysis on synthetic data, converted to SpatialExperiment
#' library(SpatialExperiment)
#' seX = SpatialExperiment(assays = list("transcripts" = t(X)), spatialCoords = Cx)
#' seY = SpatialExperiment(assays = list("metabolites" = t(Y)), spatialCoords = Ey)
#' resModtGPs = sbivar(seX, seY, assayX = "transcripts", assayY = "metabolites", method = "GPs")
#' #Multi-image analysis on Vicari data
#' data(Vicari)
#' VicariRes = sbivar(Vicari$TranscriptOutcomes, Vicari$MetaboliteOutcomes,
#' Vicari$TranscriptCoords, Vicari$MetaboliteCoords, method = "Moran", wo = "distance")
setGeneric("sbivar", function(X, ...) standardGeneric("sbivar"))
