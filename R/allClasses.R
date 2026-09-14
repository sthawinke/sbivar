setClass(
    "SbivarResults",
    slots = c(
        result = "data.frame",
        method = "character",
        multi = "logical",
        normX = "character",
        normY = "character"
    )
)
setClassUnion("numericOrNULL", c("character", "NULL"))
setClass(
    "SbivarResultsMoransI",
    contains = "SbivarResults",
    slots = c(
        maxIxy = "numericOrNULL",
        estimateSEsMoransI = "logical",
        wo = "character",
        wParams = "numeric"
    )
)
setOldClass("corSpatial")
setClassUnion("corSpatialOrNULL", c("character", "NULL"))
setClass(
    "SbivarResultsGAMs",
    contains = "SbivarResults",
    slots = c(
        families = "list",
        correlation = "corSpatialOrNULL",
        Gamm = "logical"
    )
)
