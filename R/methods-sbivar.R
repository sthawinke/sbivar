#' @rdname sbivar
#' @export
#' @param coordVars Names of coordinates
#' @param imageIdentifier A character vector of variables whose unique combinations
#' define the separate point patterns (images)
#' @param imageVars Covariates belonging to the point patterns
#' @param pointVars Names of event-wise covariates such as gene or cell for each single point
setMethod("sbivar", "MultiAssayExperiment", function(X, ...) {

})
#' @rdname sbivar
#' @export
#' @inheritParams sbivarSingle
setMethod("sbivar", "matrix", function(X, Y, Cx, Ey, ...) {
    if(!(is.matrix(Y)||is.data.frame(Y)) || !(is.matrix(Cx)||is.data.frame(Cx)) ||
       !(is.matrix(Ey)||is.data.frame(Ey))){
        stop("Since X is a matrix or dataframe, Y, Cx and Ey must be so too!")
    }
    sbivarSingle(X, Y, Cx, Ey, ...)
})
#' @rdname sbivar
#' @export
#' @importFrom spatstat.geom is.ppp
setMethod("sbivar", "list", function(X, ...) {

})
#' @param assaX,assayY Assay names to be used in the analysis
#' @rdname sbivar
#' @export
#' @importFrom SpatialExperiment SpatialExperiment spatialCoords
setMethod("sbivar", "SpatialExperiment", function(X, Y, assayX, assayY, ...) {
    if(!inherits(Y, "SpatialExperiment")){
        stop("Since X is a SpatialExperiment object, Y must be so too!")
    }
    sbivarSingle(assays(X, assayX), assays(Y, assayY), spatialCoords(X), spatialCoords(Y), ...)
})
