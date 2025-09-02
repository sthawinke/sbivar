#' Name a character vector after itself
#' @param x The vector to be names
#' @return the named vector
#' @export
#' @examples
#' selfName(LETTERS[1:5])
selfName = function(x){names(x)=x;x}
#' Convert z-value to p-value
#'
#' @param z The z-value to be converted
#' @return The p-value
makePval = function(z){
    z[is.na(z)] = 0
    tmp = pnorm(z, lower.tail = TRUE)
    tmp[z>0] = pnorm(z[z>0], lower.tail = FALSE)
    2*unname(tmp)
}
#' Scale to [0,1] range
#' @param y The vector to be scaled
#' @param na.rm passed onto \link[base]{min} and \link[base]{range}
#' @return The scaled vector
scaleZeroOne = function(y, na.rm = TRUE){
    (y-min(y, na.rm = na.rm))/diff(range(y, na.rm = na.rm))
}
#' Scale to [-1,1] range
#' @inheritParams scaleZeroOne
#' @return The scaled vector
scaleMinusOne = function(y, na.rm = TRUE){
    scaleZeroOne(y, na.rm = na.rm)*2-1
}
#' Make unique names
#' @param featX,featY vectors of feature names
#' @return A vector of names
makeNames = function(featX, featY){
    make.names(apply(expand.grid(featX, featY), 1, paste, collapse = "_"))
}
#' A wrapper for Matrix::bdiag maintaining names
#'
#' @param A,B Matrix to be used in \link[Matrix]{bdiag}
#' @return Same as \link[Matrix]{bdiag} but with dimnames
bdiagn = function(A, B){
    M <- bdiag(A, B)
    # Build new dimnames from components
    dimnames(M) <- list(c(rownames(A), rownames(B)), c(colnames(A), colnames(B)))
    M
}
#' Find trace of a matrix, of traces of an array
#'
#' A (mxm) matric has one trace (the product of the diagonal elements), a (mxmxp) array has p traces
#'
#' @param x Matrix or array
#' @param dim Dimensions defining matrices to find traces over
#'
#' @returns A trace or vector of traces
#' @importFrom methods is
tr = function(x, dim = c(1,2)) {
    if(is.matrix(x)|| is(x, "dgeMatrix") || is(x, "dgCMatrix")){
        sum(diag(x))
    } else if(is.array(x)){
        apply(x, dim, tr)
    } else {
        stop("Trace function not implemented for ", class(x))
    }
}
