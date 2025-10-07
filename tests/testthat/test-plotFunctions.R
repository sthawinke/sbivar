context("Plotting functions for single and multiple images")
test_that("Single image plotting function works", {
    expect_is(plotTopPair(resModtTest, X = X, Y = Y, Cx = Cx, Ey = Ey), "ggplot")
    expect_is(plotPairSingle( X = X, Y = Y, Cx = Cx, Ey = Ey, features = c("X1", "Y1")), "ggplot")
})
test_that("Multi-image plotting function works", {
    expect_is(plotTopPair(resGAMsMulti, Xl = Xl, Yl = Yl, Cxl = Cxl, Eyl = Eyl), "ggplot")
    expect_is(plotTopPair(resGAMsMulti, Xl, Yl, Cxl, Eyl, parameter = "cofactor"), "ggplot")
    expect_is(plotPairMulti(Xl, Yl, Cxl, Eyl, features = c("X1", "Y1")), "ggplot")
})
