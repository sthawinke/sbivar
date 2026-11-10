test_that("Array functions yield the same results as using loops", {
  m <- 3
  p <- 4
  mat <- matrix(rnorm(m^2), m, m)
  arr <- array(rnorm(m^2 * p), dim = c(m, m, p))
  expect_equal(
    arrayMatProd(arr, mat),
    vapply(seq_len(p), FUN.VALUE = mat, function(i) {
      mat %*% arr[, , i]
    })
  )
})
