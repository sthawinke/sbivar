context("Fitting single GAMs and GPs")
test_that("Fitting GAMs and GPs on single outcomes works", {
    expect_s4_class(fitGAM(data.frame("out" = X[,1], Cx), "out"), "gam")
    expect_identical(names(fitGP(X[,1], Cx, GPmethod = "REML")),
                     c("range", "nugget", "sigma", "mean"))
})
