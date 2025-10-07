context("Plot GAM fits and correlations")
test_that("GAM fitting plotting works", {
    expect_is(plotGAMs(X[, 1], Y[, 1], Cx, Ey), "ggplot")
    expect_is(plotGAMsFromMatrix(X, Y, features = c("X1", "Y1"), Cx = Cx, Ey = Ey), "ggplot")
    expect_is(plotGAMsTopResults(resGAMsSingle, X, Y, Cx = Cx, Ey = Ey), "ggplot")
})
