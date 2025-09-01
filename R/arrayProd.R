#' Multiply an array with a matrix
#' @description Array resulting from all matrix products of slices of array A mxmxp, with matrix M is mxm
arrayMatProd = function(A, M){
    n = dim(A)[1];p = dim(A)[3]
    dn = dimnames(A)
    # # Reshape A into n x (n*p)
    A <- matrix(A, n, n * p)
    # Multiply in one go
    result <- M %*% A
    # Reshape back to n x n x p
    result <- array(result, dim = c(n, n, p))
    dimnames(result) = dn
    return(result)
}
#' Multiply two arrays of same dimensions along the first two dimensions
#' @description Array resulting from all matrix products of slices of array A mxmxp
#' with those of array B of the same dimensions,
#'
arrayProd = function(A, B){
    n = dim(A)[1];p = dim(A)[3];k = dim(B)[3]
    # Reshape each into (n^2) x p
    A_mat <- matrix(A, n*n, p)
    B_mat <- matrix(B, n*n, k)
    # p x p result: matrix of all slice inner products
    C <- crossprod(A_mat, B_mat)
    dimnames(C) = list(dimnames(A)[[3]], dimnames(B)[[3]])
    C
}
#' Find traces of all inner products of matrices composing arrays A and B, yielding a pxp matrix
#'
#' @param A,B mxmxp arrays
#' @return pxp matrix of traces
arrayProd2tr <- function(A, B) {
    #High memory, but fast version
    p <- dim(A)[3]
    k <- dim(B)[3]
    n <- dim(A)[1]
    # reshape: each slice becomes a column
    Amat <- matrix(aperm(A, c(2,1,3)), n*n, p)
    Bmat <- matrix(B, n*n, k)
    # big matrix multiplication gives all traces at once
    result <- t(crossprod(Amat, Bmat)) # k x p
    dimnames(result) <- list(dimnames(B)[[3]], dimnames(A)[[3]])
    result
}
