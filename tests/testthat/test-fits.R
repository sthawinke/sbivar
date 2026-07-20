context("Fitting single GAMs and GPs")
test_that("Fitting GAMs and GPs on single outcomes works", {
    expect_s3_class(fitGAM(data.frame("out" = X[, 1], Cx), "out", Gamm = FALSE), "gam")
    expect_s3_class(fitGAM(data.frame("out" = X[, 1], Cx), "out", correlation = corExp(form = ~ x + y, nugget = TRUE, value = c(1, 0.25)), Gamm = TRUE), "gam")
    expect_named(
        fitGP(X[, 1], Cx,
            GPmethod = "REML",
            optControl = lmeControl(
                opt = "optim", maxIter = 5e2, msMaxIter = 5e2,
                niterEM = 1e3, msMaxEval = 1e3
            ),
            correlation = corGaus(form = ~ x + y, nugget = TRUE, value = c(1, 0.25))
        ),
        c("range", "nugget", "sigma", "mean")
    )
})
