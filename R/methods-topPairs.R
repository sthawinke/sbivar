#' @importFrom utils head
#' @rdname topPairs
#' @param ... passed onto \link[utils]{head}
#' @examples
#' example(sbivar, "sbivar")
#' topPairs(resGAMs)
setMethod(
    "topPairs",
    "sbivarResults",
    function(x, parameter = "Intercept", ...) {
        if (x@multi) {
            head(x@result[[parameter]], ...)
        } else {
            head(x@result, ...)
    }
})
