setClass(
    "sbivarResults",
    slots = c(
        result = "data.frame",
        method = "character",
        multi = "logical",
        normX = "character",
        normY = "character"
    )
)
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClass(
    "sbivarResultsMoransI",
    contains = "sbivarResults",
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
    "sbivarResultsGAMs",
    contains = "sbivarResults",
    slots = c(
        families = "list",
        correlation = "corSpatialOrNULL",
        Gamm = "logical"
    )
)
