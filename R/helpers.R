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
    make.names(apply(expand.grid(featX, featY), 1, paste, collapse = "__"))
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
#' A (mxm) matric has one trace (the sum of the diagonal elements), a (mxmxp) array has p traces
#'
#' @param x Matrix or array
#' @param dim Dimensions defining matrices to find traces over
#'
#' @returns A trace or vector of traces
#' @importFrom methods is
tr = function(x, dim = c(1,2)) {
    if(is.matrix(x) || is(x, "Matrix")){
        sum(diag(x))
    } else if(is.array(x)){
        apply(x, dim, tr)
    } else {
        stop("Trace function not implemented for ", class(x))
    }
}
#' Get all categrocial variables from a dataframe
#'
#' @param df The data frane
#'
#' @returns A character vector of variable names
getDiscreteVars = function(df){
    colnames(df)[!vapply(df, FUN.VALUE = TRUE, is.numeric)]
}
#' Scale variables to the [0,1] or [-1,1] range
#'
#' Scale variables to comparable ranges for plotting reasons
#' @param y The variable to scale
#' @param na.rm Should NA values be removed
#'
#' @returns scaled variable
scaleZeroOne = function(y, na.rm = TRUE){
    (y-min(y, na.rm = na.rm))/diff(range(y, na.rm = na.rm))
}
scaleMinusOne = function(y, na.rm = TRUE){
    (y-min(y, na.rm = na.rm))/diff(range(y, na.rm = na.rm))*2-1
}
#' Wrapper to normalize, select feature and scale
#'
#' @details Returns vector of NA if feature not found, leading to grey in the plots
#' @param X data matrix
#' @param feat the feature name
#' @param normFun normalizing function
#'
#' @returns A vector of values
scaleHelpFun = function(X, feat, normFun){
    normFun = match.fun(normFun)
    if(feat %in% colnames(X)){
        scaleZeroOne(normFun(X)[, feat])
    } else {
        rep(NA, nrow(X))
    }
}
#' Log normalize a data matrix
#'
#' Normalize to relative expression, add pseudocount and log-normalize
#' @param x The matrix
#' @param pseudoCount A pseudocount added to avoid taking the log of zero
#'
#' @returns A log-normalized matrix
#' @export
#' @examples
#' mat = matrix(rpois(2000), 40, 50)
#' lnMat = logNorm(mat)
logNorm = function(x, pseudoCount = 1e-8){
    stopifnot(is.matrix(x), all(x>=0))
    x = x[rowSums(x)>0,,drop = FALSE]
    dn = dimnames(x)
    out = log((as.matrix(x)+pseudoCount)/rowSums(x))
    dimnames(out) = dn
    return(out)
}
#' Split a string
#'
#' @param string The string
#' @param split string to split by
#'
#' @returns A character vector of length 2
sund = function(string, split = "__"){
    strsplit(string, split = split)[[1]]
}

#' Replace the left hand side of a formula by a fixed string
#'
#' @param x a formula
#' @param repl the replacement string
#'
#' @returns A formula
#' @importFrom stats formula
replaceLhs = function(x, repl = "out"){
    out = paste(repl, "~", deparse(x[[length(x)]]))
    return(formula(out))
}
