test_that("Helper functions do their job", {
    expect_identical(makePval(0), 1)
    expect_named(selfName("gene1"), "gene1")
    expect_true({
        foo <- range(scaleZeroOne(rnorm(20)))
        (foo[1] >= 0) && (foo[2] <= 1)
    })
    expect_true({
        foo <- range(scaleMinusOne(rnorm(20)))
        (foo[1] >= -1) && (foo[2] <= 1)
    })
    expect_identical(
        makeNames(c("gene1", "gene2"), c("compoundA", "compoundB")),
        c(
            "gene1__compoundA", "gene2__compoundA",
            "gene1__compoundB", "gene2__compoundB"
        )
    )
    mat <- matrix(rnorm(9, mean = 2), 3, 3, dimnames = list(LETTERS[1:3], letters[1:3]))
    mat2 <- matrix(rnorm(9), 3, 3, dimnames = list(
        paste0("X", LETTERS[1:3]),
        paste0("X", letters[1:3])
    ))
    expect_identical(
        dimnames(bdiagn(mat, mat2)),
        list(
            c(rownames(mat), rownames(mat2)),
            c(colnames(mat), colnames(mat2))
        )
    )
    expect_identical(as.character(replaceLhs(~ (1 | rat))[[2]]), "out")
    expect_identical(as.character(replaceLhs(outcome ~ (1 | rat))[[2]]), "out")
    expect_identical(dimnames(mat), dimnames(normMat(abs(mat), "rel")))
    expect_identical(dimnames(mat), dimnames(normMat(abs(mat), "log")))
})
