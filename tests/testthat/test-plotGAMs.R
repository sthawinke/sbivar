context("Plot GAM fits and correlations")
test_that("GAM fitting plotting works", {
    expect_is(plotGAMs(X, Y, Cx, Ey, features = c("X1", "Y2")), "ggplot")
    expect_is(plotGAMsTopResults(resGAMsSingle, X, Y, Cx = Cx, Ey = Ey), "ggplot")
    expect_is(plotGAMs(Xl, Yl, Cx, Ey, features = c("X1", "Y2")), "ggplot")
    expect_is(plotGAMsTopResults(resGAMsMulti, Xl, Yl, Cx = Cxl, Ey = Eyl), "ggplot")
})
