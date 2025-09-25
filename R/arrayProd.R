#' Multiply an array with a matrix
#'
#' @param A An mxmxp array
#' @param M An mxm matrix
#' @return An mxmxp array
#' @description Array resulting from all matrix products of slices of array A mxmxp, with matrix M is mxm.
#' This is a faster version of \code{vapply(seq_len(p), FUN.VALUE = mat, function(i)\{mat \%*\% arr[,,i]\})},
#' although it may consume more memory
#' @details The speedup comes from a single call to %*%, very efficient in BLAS.
arrayMatProd = function(A, M){
    n = dim(A)[1];p = dim(A)[3]
    dn = dimnames(A)
    # # Reshape A into n x (n*p)
    A <- matrix(A, n, n * p)
    # Multiply in one go
    C <- M %*% A
    # Reshape back to n x n x p
    C <- array(C, dim = c(n, n, p), dimnames = dn)
    return(C)
}
#' Find traces of all inner products of matrices composing arrays A and B, after transposing the second yielding a pxp matrix
#'
#' @param A,B mxmxp arrays
#' @return A pxp matrix of traces
#' @description Returns the matrix resulting from the traces of all matrix products of slices of array A mxmxp
#' with those of array B of the same dimensions. It is a faster version of
#' \code{vapply(seq_len(p), FUN.VALUE = double(p), function(i)\{
#'      vapply(seq_len(p), FUN.VALUE = double(1), function(j)\{
#'          tr(crossprod(arr2[,,i], arr[,,j]))
#'      \})
#' \})}
#'
#'
#' @details The speedup comes from a single call to crossprod, very efficient in BLAS.
arrayProdTr = function(A, B){
    n <- dim(A)[1];p = dim(A)[3];k = dim(B)[3]
    dn <- list(dimnames(A)[[3]], dimnames(B)[[3]])
    # Reshape each into (n^2) x p
    A <- matrix(A, n^2, p)
    B <- matrix(B, n^2, k)
    # p x p result: matrix of all slice inner products
    C <- crossprod(A, B)
    dimnames(C) = dn
    C
}
#' Find traces of all inner products of matrices composing arrays A and B yielding a pxp matrix
#' @description A faster version of \code{vapply(seq_len(p), FUN.VALUE = double(p), function(i)\{
#' vapply(seq_len(p), FUN.VALUE = double(1), function(j)\{
#'    tr(arr[,,i] \%*\% arr2[,,j])
#'  \})
#'\})}
#'
#' @param A,B mxmxp arrays
#' @return pxp matrix of traces
arrayProd2tr <- function(A, B) {
    #High memory, but fast version
    p <- dim(A)[3]
    k <- dim(B)[3]
    n <- dim(A)[1]
    dn <- list(dimnames(B)[[3]], dimnames(A)[[3]])
    # reshape: each slice becomes a column
    A <- matrix(aperm(A, c(2,1,3)), n*n, p)
    B <- matrix(B, n*n, k)
    # big matrix multiplication gives all traces at once
    C <- t(crossprod(A, B)) # k x p
    dimnames(C) <- dn
    C
}
