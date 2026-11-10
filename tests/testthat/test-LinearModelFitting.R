test_that("fitLinModels works for GAM input", {
  expect_s4_class(
    type = "sbivarResults",
    multiFitGams <- fitLinModels(estGAMs,
      designDf = toyDesign, inverseWeigh = TRUE,
      Formula = out ~ covariate + cofactor + (1 | group)
    )
  )
    expect_s4_class(
        type = "sbivarResults",
    multiFitGams2 <- fitLinModels(estGAMs,
      inverseWeigh = FALSE,
      designDf = toyDesign,
      Formula = out ~ covariate + cofactor + (1 | group)
    )
  )
    expect_s4_class(
        type = "sbivarResults",
    multiFitGams3 <- fitLinModels(estGAMs,
      inverseWeigh = FALSE,
      Formula = out ~ 1
    )
  )
  # Extract the results
  expect_s4_class(resGams <- extractResultsMulti(multiFitGams, designDf = toyDesign), "sbivarResults")
  expect_named(
    resGams@result,
    c("Intercept", "covariate", "cofactor")
  )
  expect_identical(colnames(resGams@result$Intercept), c("Modality_X", "Modality_Y", "Estimate", "SE", "pVal", "pAdj"))
  expect_warning(
    resGams3 <- extractResultsMulti(multiFitGams3)
  )
  expect_named(resGams3@result, "Intercept")
  expect_false(is.unsorted(resGams@result$Intercept[, "pVal"]))
})
test_that("fitLinModels works for Moran's I input", {
    expect_s4_class(
        type = "sbivarResults",
    multiFitMoran <- fitLinModels(estMoran,
      designDf = toyDesign,
      Formula = out ~ covariate + cofactor + (1 | group)
    )
  )
    expect_s4_class(
        type = "sbivarResults",
    multiFitMoran2 <- fitLinModels(estMoran,
      inverseWeigh = FALSE,
      designDf = toyDesign,
      Formula = out ~ covariate + cofactor + (1 | group)
    )
  )
  expect_type(
    type = "list",
    multiFitMoran3 <- fitLinModels(estMoran,
      inverseWeigh = FALSE, scaleByMax = FALSE,
      designDf = toyDesign,
      Formula = out ~ covariate + cofactor + (1 | group)
    )
  )
  expect_s4_class(
      type = "sbivarResults",
    multiFitMoran4 <- fitLinModels(estMoran,
      inverseWeigh = TRUE, scaleByMax = FALSE,
      designDf = toyDesign,
      Formula = out ~ covariate + cofactor + (1 | group)
    )
  )
})
test_that("fitLinModels works for correlation input", {
    expect_s4_class(
        type = "sbivarResults",
    multiFitCors <- fitLinModels(estMultiCor,
      designDf = toyDesign,
      Formula = out ~ covariate + cofactor + (1 | group)
    )
  )
})
test_that("Model fitting and extraction throws errors where appropriate", {
  expect_error(fitLinModels(estMoran, Formula = out ~ covariate + cofactor + (1 | group)))
})
