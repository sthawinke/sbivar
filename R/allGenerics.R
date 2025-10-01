#' Perform a spatial bivariate association analysis
#' @description
#' Perform the spatial analysis of choice, either on a single or multiple imgaes
#' @rdname sbivar
#' @importFrom methods setGeneric setMethod
#' @export
#' @param X the input object, see methods('sbivar')
#' @param ... additional constructor arguments
#' @return The analysis result
#' @seealso \link[spatstat.geom]{hyperframe}
#' @examples
#' data(Vicari)
#' sbivarSingle(singleStx, singleMet, singleStxCoords, singleMetCoords,
#' method = "GAMs",  families = FamiliesViari)
setGeneric("sbivar", function(X, ...) standardGeneric("sbivar"))
