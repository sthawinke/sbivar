context("Unit tests for input errors")
test_that("SbivarSingle works for correct input", {
    expect_is(sbivarSingle(X, Y, Cx, Ey, method = "GAMs"), "matrix")
    expect_identical(colnames(sbivarSingle(X, Y, Cx, Ey, method = "Modified")),
                     c("Correlation", "Effective sample size", "pVal", "pAdj"))
    expect_false(is.unsorted(sbivarSingle(X, Y[seq_len(nrow(X)),], Cx,
                                          method = "Modified")[, "pVal"]))
    expect_message(sbivarSingle(X, X, Cx, Cx, method = "Modified"))
    expect_silent(sbivarSingle(X, Y, Cx, Ey, method = "GAMs"))
    expect_is(sbivarSingle(X, Y, Cx, Ey, method = "GPs"), "matrix")
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
    expect_error(sbivarSingle(X, Y, X, Ey, method = "Modified"))
    expect_error(sbivarSingle(X, Y, Cx[-1,], Ey, method = "GAMs"))
    expect_error(sbivarSingle(X, Y, Cx, Ey, method = "GPs",
        corStruct = nlme::corExp(form = ~ x + y, nugget = TRUE, value = c(1, 0.25))))

})
test_that("sbivarMulti works for correct input", {
    gamList <- sbivarMulti(Xl, Yl, Cxl, Eyl, method = "GAMs")
    expect_is(gamList, "list")
    expect_identical(names(gamList), c("estimates", "method"))
    expect_is(gamList$estimates[[1]], "matrix")
    expect_is(sbivarMulti(Xl, Yl, Cxl, Eyl, method = "Moran's I"), "list")
    expect_is(sbivarMulti(Xl, Yl, Cxl, Eyl, method = "Correlation"), "list")
})
test_that("sbivarMulti throws errors for incorrect input", {
    expect_error(sbivarMulti(Xl, Yl, Cxl, Eyl, method = "GPs"))
    expect_error(sbivarMulti(Xl, Yl, Cxl, Eyl, method = "Modified t-test"))
    expect_error(sbivarMulti(Xl, Yl, Cxl, method = "GAMs"))
    expect_error(sbivarMulti(Xl, Yl, Cxl, method = "GAMs",
                             families = list("X" = gaussian(), "Y" = mean)))
    expect_error(sbivarMulti(Xl, Yl, Cxl[-1], Eyl, method = "Moran's I"))
    expect_error(sbivarMulti(Xl, Yl, Cxl, Cxl, method = "Moran's I"))
})
