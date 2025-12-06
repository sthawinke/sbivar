context("Cauychy combination rule")
test_that("Cauchy combination rule works as expected", {
    pVals = runif(20, 0 ,1)
    cctP = CCT(pVals)
    expect_is(cctP, "numeric")
    expect_true((cctP <= 1) && (cctP >= 0))
    expect_warning(CCT(c(1, pVals)))
    expect_is(c(cctP, 1e-20), "numeric")
})

