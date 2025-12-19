context("Plot coordinate matrices")
test_that("Coordinate plotting works", {
    expect_silent(plotCoords(Cx, Ey))
    expect_silent(plotCoordsMulti(Cxl, Eyl))
})
