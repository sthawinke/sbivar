context("Unit tests for array products and traces")
test_that("Array functions yield the same results as using loops", {
    m = 3;p = 4
    mat = matrix(rnorm(m^2), m, m)
    arr = array(rnorm(m^2*p), dim = c(m,m,p))
    mat2 = matrix(rnorm(m^2), m, m)
    arr2 = array(rnorm(m^2*p), dim = c(m,m,p))
    expect_equal(arrayMatProd(arr, mat),
                     vapply(seq_len(p), FUN.VALUE = mat, function(i){
                          mat %*% arr[,,i]
                     }))
    expect_equal(unname(arrayProdTr(arr, arr2)),
                     vapply(seq_len(p), FUN.VALUE = double(p), function(i){
                         vapply(seq_len(p), FUN.VALUE = double(1), function(j){
                            tr(crossprod(arr2[,,i], arr[,,j]))
                         })
                     }))
    expect_equal(unname(arrayProd2tr(arr, arr2)),
                 vapply(seq_len(p), FUN.VALUE = double(p), function(i){
                     vapply(seq_len(p), FUN.VALUE = double(1), function(j){
                         tr(arr[,,i] %*% arr2[,,j])
                     })
                 }))
})
