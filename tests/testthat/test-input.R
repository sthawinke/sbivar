context("Unit tests for input errors")
test_that("SbivarSingle works for correct input", {
    expect_is(sbivar(X, Y, Cx, Ey, method = "GAMs"), "matrix")
    expect_identical(colnames(sbivar(X, Y, Cx, Ey, method = "Modified")),
                     c("Correlation", "Effective sample size", "pVal", "pAdj"))
    expect_false(is.unsorted(sbivar(X, Y[seq_len(nrow(X)),], Cx,
                                          method = "Modified")[, "pVal"]))
    expect_message(sbivar(X, X, Cx, Cx, method = "Modified"))
    expect_silent(sbivar(X, Y, Cx, Ey, method = "GAMs", verbose = FALSE))
    expect_is(sbivar(X, Y, Cx, Ey, method = "GPs"), "matrix")
})
test_that("SbivarSingle throws errors for incorrect input", {
    expect_error(sbivar(X, Y, Cx, Ey, method = "Moran's I"))
    expect_error(sbivar(X, Y, Cx, method = "GAMs"))
    expect_error(sbivar(X, Y, Cx, Cx, method = "GAMs"))
    expect_error(sbivar(X, Y, Cx, Ey, method = "GAMs",
                              families = list(gaussian(), gaussian())))
    expect_error(sbivar(X, Y, Cx, Ey, method = "GAMs",
                              families = list("X" = gaussian(), "Y" = mean)))
    expect_error(sbivar(X, Y, Cx, Y, method = "Modified"))
    expect_error(sbivar(X, Y, X, Ey, method = "Modified"))
    expect_error(sbivar(X, Y, Cx[-1,], Ey, method = "GAMs"))
    expect_error(sbivar(X, Y, Cx, Ey, method = "GPs",
        corStruct = nlme::corExp(form = ~ x + y, nugget = TRUE, value = c(1, 0.25))))

})
test_that("sbivarMulti works for correct input", {
    gamList <- sbivar(Xl, Yl, Cxl, Eyl, method = "GAMs")
    expect_is(gamList, "list")
    expect_identical(names(gamList), c("estimates", "method"))
    expect_is(gamList$estimates[[1]], "matrix")
    expect_is(sbivar(Xl, Yl, Cxl, Eyl, method = "Moran's I"), "list")
    expect_is(sbivar(Xl, Yl, Cxl, Eyl, method = "Correlation"), "list")
})
test_that("sbivarMulti throws errors for incorrect input", {
    expect_error(sbivar(Xl, Yl, Cxl, Eyl, method = "GPs"))
    expect_error(sbivar(Xl, Yl, Cxl, Eyl, method = "Modified t-test"))
    expect_error(sbivar(Xl, Yl, Cxl, method = "GAMs"))
    expect_error(sbivar(Xl, Yl, Cxl, method = "GAMs",
                             families = list("X" = gaussian(), "Y" = mean)))
    expect_error(sbivar(Xl, Yl, Cxl[-1], Eyl, method = "Moran's I"))
    expect_error(sbivar(Xl, Yl, Cxl, Cxl, method = "Moran's I"))
})
test_that("Sbivar works on BioConductor objects SpatialExperiment and MultiAssayExperiment", {
    library(SpatialExperiment)
    library(MultiAssayExperiment)
    n=8e1;m=1e2;p=4;k=3
    rna_counts <- matrix(rpois(n*p, lambda = 5), nrow = p, ncol = n)
    colnames(rna_counts) <- paste0("rnaSpot", seq_len(n))
    rna_coords <- cbind("x" = runif(n, 0, 1), "y" = runif(n, 0, 1))
    spe_rna <- SpatialExperiment(assays = list("counts" = rna_counts),
                                 spatialCoords = rna_coords)
    # --- Protein data ---
    prot_counts <- matrix(rpois(m*k, lambda = 10), nrow = k, ncol = m)
    colnames(prot_counts) <- paste0("protSpot", seq_len(m))
    prot_coords <- cbind("x" = runif(m, 0, 1), "y" = runif(m, 0, 1))
    spe_prot <- SpatialExperiment(assays = list("counts" = prot_counts),
                                  spatialCoords = prot_coords)
    # --- Combine into MultiAssayExperiment ---
    mae <- MultiAssayExperiment(experiments = list(RNA = spe_rna, Protein = spe_prot))
    expect_is(sbivar(mae, experiments = c("RNA", "Protein"), assayX = "counts",
             assayY = "counts", families = list("X" = mgcv::nb(), "Y" = mgcv::nb())), "matrix")
    #Provide separate SpatialExperiment objects
    expect_is(sbivar(mae[["RNA"]], mae[["Protein"]], assayX = "counts",
            assayY = "counts", families = list("X" = mgcv::nb(), "Y" = mgcv::nb())), "matrix")
})
