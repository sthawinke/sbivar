setMethod(
    "topPairs",
    "sbivarResults",
    function(x, parameter, ...) {
        if(x@multi){
            head(x@result[[parameter]], ...)
        } else {
            head(x@result, ...)
        }
})
