setClass(
    "SbivarResults",
    slots = c(
        result = "matrix",
        method = "character",
        multi = "logical",
        normX = "characeter",
        normY = "character",
    ),
    prototype = list(
        maxIxy = NA,
        wo = "",
        wParams = NA,


    )
)
setClass(
    "SbivarResultsGAM",
    contains = "SbivarResults",
    slots = c(
        families = character(2),
        correlation = FALSE,

    )
)
