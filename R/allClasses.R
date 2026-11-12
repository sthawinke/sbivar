setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("listOrDataFrameOrNull", c("data.frame", "list", "NULL"))
#' sbivarResults class
#'
#' An S4 class representing the result of a sbivar call, testing for bivariate
#' spatial association in one or more images
#'
#' @slot method The statistical test used
#' @slot normX Normalization applied to the first modality
#' @slot normY Normalization applied to the second modality
#' @slot multi A boolean, were multiple images tested?
#' @slot assayX Name of the assay of the first modality
#' @slot assayY Name of the assay of the second modality
#' @slot models Optional, the linear model fits if multi=TRUE
#' @slot result The result of the statistical test: a dataframe
#' @slot estimates Raw estimates of association measures,
#' for input to linear model fitting
#'
#' @export
setClass(
    "sbivarResults",
    slots = c(
        result = "listOrDataFrameOrNull",
        method = "character",
        normX = "character",
        normY = "character",
        multi = "logical",
        assayX = "characterOrNULL",
        assayY = "characterOrNULL",
        models = "listOrNULL",
        estimates = "listOrNULL"
    )
)
setClassUnion("numericOrMatrixOrNULL", c("numeric", "matrix", "NULL"))
setClass(
    "sbivarResultsMoransI",
    contains = "sbivarResults",
    slots = c(
        maxIxy = "numericOrMatrixOrNULL",
        returnSEsMoransI = "logical",
        wo = "character",
        wParams = "numeric"
    )
)
methods::setOldClass(c("corExp", "corSpatial", "corStruct"))
methods::setOldClass(c("corGaus", "corSpatial", "corStruct"))
methods::setOldClass(c("corLin", "corSpatial", "corStruct"))
methods::setOldClass(c("corRatio", "corSpatial", "corStruct"))
methods::setOldClass(c("corSpher", "corSpatial", "corStruct"))
setClassUnion("corSpatialOrNULL", c("corSpatial", "NULL"))
setClass(
    "sbivarResultsGAMs",
    contains = "sbivarResults",
    slots = c(
        families = "list",
        correlation = "corSpatialOrNULL",
        Gamm = "logical"
    )
)
setValidity("sbivarResults", function(object) {
    if (length(object@normX) != 1 || !is.character(object@normX)) {
        return("normX must have length 1 and be a character")
    }
    if (length(object@normY) != 1 || !is.character(object@normY)) {
        return("normY must have length 1 and be a character")
    }
    if (length(object@method) != 1 ||
        !is.character(object@method)) {
        return("method must have length 1 and be a character")
    }
    if (!is.null(object@assayX) && !(length(object@assayX) == 1 &&
        is.character(object@assayX))) {
        return("assayX must be null or have length 1 and be a character")
    }
    if (!is.null(object@assayY) && !(length(object@assayY) == 1 &&
        is.character(object@assayY))) {
        return("assayY must be null or have length 1 and be a character")
    }
    if (!is.logical(object@multi)) {
        return("multi must be a boolean")
    }
    if (!(
        is.data.frame(object@result) || is.list(object@result) ||
            is.null(object@result)
    )) {
        return("result must be a dataframe, list or null")
    }
    if (!is.data.frame(object@result) && !object@multi) {
        return("result must be a dataframe for single image analysis")
    }
    if (!(is.list(object@estimates) || is.null(object@estimates))) {
        return("estimates must be a list or null")
    }
    if (!object@multi &&
        (is.list(object@models) || is.list(object@estimates))) {
        return("Estimates or models can only be attached for multi-image analysis")
    }
    TRUE
})
#' Show basic properties of a sbivarResults object
#' @importFrom methods show
#' @param object The sbivarResults object
#' @return Prints a description to console
#' @export
setMethod("show", signature(object = "sbivarResults"), function(object) {
    cat("Sbivar fit with method", object@method, "on", if (object@multi) {
        c(length(object@result), "images")
    } else {
        c(
            "single image with",
            nrow(object@result),
            "feature pairs tested"
        )
    })
})
