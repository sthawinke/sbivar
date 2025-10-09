context("Fit linear models and extract results")
test_that("SbivarSingle works for correct input", {
    expect_is(class = "list",
              multiFitGams <- fitLinModels(estGAMs, designDf = toyDesign,
                      Formula = out ~ covariate + cofactor + (1|group)))
    expect_is(class = "list",
              multiFitMoran <- fitLinModels(estMoran, designDf = toyDesign,
                      Formula = out ~ covariate + cofactor + (1|group)))
    expect_is(class = "list",
              multiFit <- fitLinModels(estCorrelations, designDf = toyDesign,
                      Formula = out ~ covariate + cofactor + (1|group)))
    #Extract the results
    expect_identical(names(resGams <- extractResultsMulti(multiFitGams, designDf = toyDesign)$result),
                     c("Intercept", "covariate", "cofactor"))
    expect_identical(colnames(resGams$Intercept), c("Estimate", "SE", "pVal", "pAdj"))
    expect_false(is.unsorted(resGams$Intercept[, "pVal"]))
})
test_that("Model fitting and extraction throws errors where appropriate", {
    expect_error(fitLinModels(estMoran, Formula = out ~ covariate + cofactor + (1|group)))
})


