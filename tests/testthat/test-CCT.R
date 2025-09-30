context("Cauychy combination rule")
test_that("Cauchy combination rule works as expected", {
    pVals = runif(20, 0 ,1)
    cctP = CCT(pVals)
    expect_is(cctP, "numeric")
    expect_true((cctP <= 1) && (cctP >= 0))
})
test_that("Cauchy combination rule throws errors where needed", {
    pVals = runif(20, 1 ,2)
    expect_error(CCT(pVals))
})
