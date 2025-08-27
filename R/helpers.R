#' Name a character vector after itself
selfName = function(x){names(x)=x;x}
#' Convert z-value to p-value
makePval = function(z){
    z[is.na(z)] = 0
    tmp = pnorm(z, lower.tail = TRUE)
    tmp[z>0] = pnorm(z[z>0], lower.tail = FALSE)
    2*unname(tmp)
}
