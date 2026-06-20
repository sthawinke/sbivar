context("Plot spot sums")
test_that("Coordinate plotting works", {
    expect_silent(plotSpotSums(Xl, Yl, Cxl, Eyl))
})
