#' @importFrom utils head
#' @rdname topPairs
#' @param ... passed onto \link[utils]{head}
#' @param parameter For multi-image analysis, the parameter for which top
#'   results should be shown
#' @param digits The number of significant digits to be rounded to
#' @examples
#' example(sbivar, "sbivar")
#' topPairs(resGAMs)
setMethod(
    "topPairs",
    "sbivarResults",
    function(x, parameter = "Intercept", digits = 5L, ...) {
        out <- if (x@multi) {
            head(x@result[[parameter]], ...)
        } else {
            head(x@result, ...)
        }
        for(i in setdiff(names(out), c("Modality_X", "Modality_Y"))){
            out[[i]] <- signif(out[[i]], digits = digits)
        }
        print(out, row.names = FALSE)
    }
)
