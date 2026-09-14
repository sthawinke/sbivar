setClass(
    "SbivarResults",
    slots = c(
        result = "matrix",
        method = "character",
        multi = "logical",
        normX = "character",
        normY = "character"
    )
)
setClass(
    "SbivarResultsMoransI",
    contains = "SbivarResults",
    slots = c(
        maxIxy = "numeric",
        estimateSEsMoransI = "logical",
        wo = "character",
        wParams = "numeric"
    )
)
setClass(
    "SbivarResultsGAM",
    contains = "SbivarResults",
    slots = c(
        families = "character",
        correlation = "logical"
    )
)
