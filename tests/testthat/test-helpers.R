context("Unit tests for small auxiliary function")
test_that("sund yields correct result", {
    expect_identical(sund("gene1--gene2"), c("gene1", "gene2"))
    expect_identical("gene1", "gene1")
})
