setMethod(
    "plotTopPair",
    "sbivarResults",
    function(x, topRank = 1, parameter = "Intercept", scaleBySampleSums = FALSE, ...) {
        stopifnot(is.numeric(topRank), topRank >= 1, is.logical(scaleBySampleSums), is.character(parameter))
        if (x@multi) {
            stopifnot(parameter %in% names(result))
            topFeats <- x@results[[parameter]][topRank, c("Modality_X", "Modality_Y")]
            plotPairMulti(
                features = topFeats, assayX = results@assayX,
                assayY = results@assayY, normX = results@normX, scaleBySampleSums = scaleBySampleSums,
                normY = results@normY, ...
            )
        } else {
            topFeats <- x@results[topRank, c("Modality_X", "Modality_Y")]
            plotPairSingle(
                features = topFeats, assayX = results@assayX, scaleBySampleSums = scaleBySampleSums,
                assayY = results@assayY, normX = results@normX, normY = results@normY, ...
            )
        }
    }
)
