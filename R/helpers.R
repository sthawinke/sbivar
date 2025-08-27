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
