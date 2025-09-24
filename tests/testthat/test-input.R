context("Unit tests for input errors")
test_that("SbivarSingle works for correct input", {
    expect_is(sbivarSingle(X, Y, Cx, Ey, method = "GAMs"), "matrix")
    expect_identical(colnames(sbivarSingle(X, Y, Cx, Ey, method = "Modified")),
                     c("Correlation", "Effective sample size", "pVal", "pAdj"))
    expect_false(is.unsorted(sbivarSingle(X, Y[seq_len(nrow(X)),], Cx,
                                          method = "Modified")[, "pVal"]))
    expect_message(sbivarSingle(X, X, Cx, Cx, method = "Modified"))
    expect_silent(sbivarSingle(X, Y, Cx, Ey, method = "GAMs"))
})
test_that("SbivarSingle throws errors for incorrect input", {
    expect_error(sbivarSingle(X, Y, Cx, Ey, method = "Moran's I"))
    expect_error(sbivarSingle(X, Y, Cx, method = "GAMs"))
    expect_error(sbivarSingle(X, Y, Cx, Cx, method = "GAMs"))
    expect_error(sbivarSingle(X, Y, Cx, Ey, method = "GAMs",
                              families = list(gaussian(), gaussian())))
    expect_error(sbivarSingle(X, Y, Cx, Ey, method = "GAMs",
                              families = list("X" = gaussian(), "Y" = mean)))
    expect_error(sbivarSingle(X, Y, Cx, Y, method = "Modified"))
})
