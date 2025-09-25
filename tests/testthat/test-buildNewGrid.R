context("Build new grid for GAMS")
gridSize <- 1e2
test_that("buildNewGrid has correct dimensions", {
 expect_identical(dim(buildNewGrid(Cx, Ey, n_points_grid = gridSize)), c(gridSize,2))
})
