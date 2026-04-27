context("Weights exploration")
test_that("Weight plotting function works", {
    expect_silent(exploreWeights(10^c(-5, -4, -3)))
})
test_that("Weight plotting function throws errors where appropriate", {
    expect_error(exploreWeights(10^c(-5, -4, -3)), dists = 0:10)
})
