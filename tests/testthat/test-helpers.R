context("Unit tests for small auxiliary function")
test_that("Unit tests work", {
    expect_identical(makePval(0), 1)
    expect_identical(names(selfName("gene1")), "gene1")
    expect_true({
        foo = range(scaleZeroOne(rnorm(20)))
        (foo[1] >= 0) && (foo[2] <= 1)
    })
    expect_true({
        foo = range(scaleMinusOne(rnorm(20)))
        (foo[1] >= -1) && (foo[2] <= 1)
    })
    expect_identical(makeNames(c("gene1", "gene2"), c("compoundA", "compoundB")),
                     c("gene1_compoundA", "gene2_compoundA", "gene1_compoundB", "gene2_compoundB"))
})
