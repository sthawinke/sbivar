setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("dataframeOrList", c("data.frame", "list"))
setClass(
    "sbivarResults",
    slots = c(
        result = "dataframeOrList",
        method = "character",
        multi = "logical",
        normX = "character",
        normY = "character",
        assayX = "characterOrNULL",
        assayY = "characterOrNULL"
    )
)
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClass(
    "sbivarResultsMoransI",
    contains = "sbivarResults",
    slots = c(
        maxIxy = "numericOrNULL",
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
