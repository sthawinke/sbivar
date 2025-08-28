#' Name a character vector after itself
#' @param x The vector to be names
selfName = function(x){names(x)=x;x}
#' Convert z-value to p-value
#'
#' @param z The z-value to be converted
makePval = function(z){
    z[is.na(z)] = 0
    tmp = pnorm(z, lower.tail = TRUE)
    tmp[z>0] = pnorm(z[z>0], lower.tail = FALSE)
    2*unname(tmp)
}
#' Scale to [0,1] range
#' @param y The vector to be scaled
#' @param na.rm passed onto \link[stats]{min} and \link[stats]{range}
scaleZeroOne = function(y, na.rm = TRUE){
    (y-min(y, na.rm = na.rm))/diff(range(y, na.rm = na.rm))
}
#' Scale to [-1,1] range
#' @inheritParams scaleZeroOne
scaleMinusOne = function(y, na.rm = TRUE){
    scaleZeroOne(y, na.rm = na.rm)*2-1
}
#' Make unique names
#' @param featX,featY vectors of feature names
makeNames = function(featX, featY){
    make.names(apply(expand.grid(featX, featY), 1, paste, collapse = "_"))
}
