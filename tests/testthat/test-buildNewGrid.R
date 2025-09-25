context("Build new grid for GAMS")
gridSize <- 1e2
test_that("buildNewGrid has correct dimensions", {
    expect_identical(colnames(buildNewGrid(Cx, Ey, n_points_grid = gridSize)),
                     c("x", "y"))
})
