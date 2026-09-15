#' @rdname plotTopPair
setMethod(
    "plotTopPair",
    "sbivarResults",
    function(x, ..., topRank, parameter, scaleBySampleSums, normX, normY) {
        stopifnot(is.numeric(topRank), topRank >= 1, is.logical(scaleBySampleSums), is.character(parameter))
        if (x@multi) {
            stopifnot(parameter %in% names(x@result))
            topFeats <- x@result[[parameter]][topRank, c("Modality_X", "Modality_Y")]
            plotPairMulti(
                features = topFeats, assayX = x@assayX,
                assayY = x@assayY, normX = normX, scaleBySampleSums = scaleBySampleSums,
                normY = normY, ...
            )
        } else {
            topFeats <- x@result[topRank, c("Modality_X", "Modality_Y")]
            plotPairSingle(
                features = topFeats, assayX = x@assayX, scaleBySampleSums = scaleBySampleSums,
                assayY = x@assayY, normX = normX, normY = normY, ...
            )
        }
    }
)
