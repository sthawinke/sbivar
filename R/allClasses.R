setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("dataframeOrList", c("data.frame", "list"))
setClass(
  "sbivarResults",
  slots = c(
    result = "dataframeOrList",
    method = "character",
    normX = "character",
    normY = "character",
    multi = "logical",
    assayX = "characterOrNULL",
    assayY = "characterOrNULL",
    models = "listOrNULL"
  )
)
setClassUnion("numericOrmatrixOrNULL", c("numeric", "matrix", "NULL"))
setClass(
  "sbivarResultsMoransI",
  contains = "sbivarResults",
  slots = c(
    maxIxy = "numericOrmatrixOrNULL",
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
#' @export
setMethod(
  "show",
  signature(object = "sbivarResults"),
  function(object) {
    cat(
      "Sbivar fit with ", object@method, " method on",
      if (object@multi) {
        c(length(object@result), "images")
      } else {
        c("single image with", nrow(object@result), "feature pairs tested")
      }
    )
  }
)
setValidity("sbivarResults", function(object) {
  if (length(object@normX) != 1 || !is.character(object@normX)) {
    return("normX must have length 1 and be a character")
  }
  if (length(object@normY) != 1 || !is.character(object@normY)) {
    return("normY must have length 1 and be a character")
  }
  if (length(object@method) != 1 || !is.character(object@method)) {
    return("method must have length 1 and be a character")
  }
  if (!is.null(object@assayX) && !(length(object@assayX) == 1 && is.character(object@assayX))) {
    return("assayX must be null or have length 1 and be a character")
  }
  if (!is.null(object@assayY) && !(length(object@assayY) == 1 && is.character(object@assayY))) {
    return("assayY must be null or have length 1 and be a character")
  }
  if (!is.logical(object@multi)) {
    return("multi must be a boolean")
  }
  if (!object@multi && is.list(object@models)) {
        return("Models can only be attached for multi-image analysis")
   }
  TRUE
})
