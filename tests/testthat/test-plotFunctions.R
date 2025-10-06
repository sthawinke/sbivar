context("Plotting functions for single and multiple images")
test_that("Single image plotting function works", {
    expect_is(plotTopPairSingle(resModtTest, X, Y, Cx, Ey),
              "ggplot")
    expect_is(plotPairSingle(X, Y, Cx, Ey, features = c("X1", "Y1")), "ggplot")
})
test_that("Multi-image plotting function works", {
    expect_is(plotTopPairMulti(resGams, Xl, Yl, Cxl, Eyl), "ggplot")
    expect_is(plotTopPairMulti(resGams, Xl, Yl, Cxl, Eyl, parameter = "cofactor"), "ggplot")
    expect_is(plotPairMulti(Xl, Yl, Cxl, Eyl, features = c("X1", "Y1")), "ggplot")
})
