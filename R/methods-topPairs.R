setMethod(
    "topPairs",
    "sbivarResults",
    function(x, ...) head(x@result, ...)
)
