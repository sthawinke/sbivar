context("Unit tests for input errors")
test_that("SbivarSingle works for correct input", {
    sbiRes <- sbivar(X, Y, Cx, Ey, method = "GAMs")
    expect_is(sbiRes, "list")
    sbiResNoGp <- sbivar(X, Y, Cx, Ey, method = "GAMs", includeGPsmooth = FALSE)
    expect_is(sbiResNoGp, "list")
    sbiResGPtest <- sbivar(X, Y, Cx, Ey, method = "GAMs", testSmooth = "field")
    expect_is(sbiResGPtest, "list")
    expect_message(sbiResMoran <- sbivar(X, Y, Cx, Ey, method = "Moran's I"))
    expect_message(sbiResMod <- sbivar(X, X, Cx, method = "Modified"))
    expect_true(all(c("pVal", "pAdj") %in% colnames(sbiResMoran$result)))
    expect_false(is.unsorted(sbiRes$result[, "pVal"]))
    expect_silent(sbivar(X, Y, Cx, Ey, method = "GAMs", verbose = FALSE))
    expect_is(sbivar(X, Y, Cx, Ey, method = "GPs"), "list")
})
test_that("SbivarSingle throws errors for incorrect input", {
    expect_error(sbivar(X, Y, cbind(Cx, Cx), Ey, method = "GAMs"))
    expect_error(sbivar(X, Y, Cx, Ey, method = "Modified"))
    expect_error(sbivar(X, Y, Cx + 10, Ey, method = "GAMs"))
    expect_error(sbivar(X, Y, Cx, Ey + 10, method = "GAMs"))
    expect_error(sbivar(X, Y, Cx, method = "GAMs"))
    expect_error(sbivar(X, Y, Cx, Cx, method = "GAMs"))
    expect_error(sbivar(X, Y, Cx, Ey, method = "GAMs", testSmooth = "all"))
    expect_error(sbivar(X, Y, Cx, Ey,
        method = "GAMs",
        families = list(gaussian(), gaussian())
    ))
    expect_error(sbivar(X, Y, Cx, Ey,
        method = "GAMs",
        families = list("X" = gaussian(), "Y" = mean)
    ))
    expect_error(sbivar(X, Y, Cx[-1, ], Ey, method = "GAMs"))
    expect_error(sbivar(X, Y, Cx, Ey,
        method = "GPs",
        corStruct = nlme::corExp(form = ~ x + y, nugget = TRUE, value = c(1, 0.25))
    ))
    Xunnamed <- X
    colnames(Xunnamed) <- NULL
    expect_error(sbivar(Xunnamed, Y, Cx, Ey, method = "GAMs"))
    Xunnamed <- X
    rownames(Xunnamed) <- NULL
    expect_error(sbivar(Xunnamed, Y, Cx, Ey, method = "GAMs"))
})
test_that("sbivarMulti works for correct input", {
    gamList <- sbivar(Xl, Yl, Cxl, Eyl, method = "GAMs")
    expect_is(gamList, "list")
    expect_identical(names(gamList), c("estimates", "method", "multi", "normX", "normY", "families"))
    expect_is(gamList$estimates[[1]]$res, "matrix")
    expect_is(sbivar(Xl, Yl, Cxl, Eyl, method = "Moran's I"), "list")
})
test_that("sbivarMulti throws errors for incorrect input", {
    expect_error(sbivar(Xl, Yl, Cxl, Eyl, method = "GPs"))
    expect_error(sbivar(Xl, Yl, Cxl, Eyl, method = "Modified t-test"))
    expect_error(sbivar(Xl, Yl, Cxl, method = "GAMs"))
    expect_error(sbivar(Xl, Yl, Cxl,
        method = "GAMs",
        families = list("X" = gaussian(), "Y" = mean)
    ))
    expect_error(sbivar(Xl, Yl, Cxl[-1], Eyl, method = "Moran's I"))
    expect_error(sbivar(Xl, Yl, Cxl, Cxl, method = "Moran's I"))
})
library(SpatialExperiment)
library(MultiAssayExperiment)
n <- 8e1
m <- 1e2
p <- 4
k <- 3
rna_counts <- matrix(rpois(n * p, lambda = 5),
                     nrow = p, ncol = n,
                     dimnames = list(paste0("rnaSpot", seq_len(p)), paste0("Gene", seq_len(n)))
)
rna_coords <- cbind("x" = runif(n, 0, 1), "y" = runif(n, 0, 1))
rownames(rna_coords) <- colnames(rna_counts)
spe_rna <- SpatialExperiment(
    assays = list("counts" = rna_counts),
    spatialCoords = rna_coords
)
# --- Protein data ---
prot_counts <- matrix(rpois(m * k, lambda = 10),
                      nrow = k, ncol = m,
                      dimnames = list(paste0("protSpot", seq_len(k)), paste0("Protein", seq_len(m)))
)
prot_coords <- cbind("x" = runif(m, 0, 1), "y" = runif(m, 0, 1))
rownames(prot_coords) <- colnames(prot_counts)
spe_prot <- SpatialExperiment(
    assays = list("counts" = prot_counts),
    spatialCoords = prot_coords
)
# --- Combine into MultiAssayExperiment ---
mae <- MultiAssayExperiment(experiments = list(RNA = spe_rna, Protein = spe_prot))
test_that("Sbivar works on BioConductor objects SpatialExperiment and MultiAssayExperiment", {
    expect_is(sbivar(mae,
        experimentX = "RNA", experimentY = "Protein", assayX = "counts",
        assayY = "counts", families = list("X" = mgcv::nb(), "Y" = mgcv::nb())
    ), "list")
    # Provide separate SpatialExperiment objects
    expect_is(sbivar(mae[["RNA"]], mae[["Protein"]],
        assayX = "counts",
        assayY = "counts", families = list("X" = mgcv::nb(), "Y" = mgcv::nb())
    ), "list")
})
test_that("Sbivar fails on BioConductor objects SpatialExperiment and MultiAssayExperiment where appropriate", {
    expect_error(sbivar(mae,
                     experimentX = "RNA", experimentY = "Protein", assayX = "counts",
                     families = list("X" = mgcv::nb(), "Y" = mgcv::nb())
    ))
    # Provide separate SpatialExperiment objects
    expect_error(sbivar(mae, mae[["Protein"]],
                     assayX = "counts",
                     assayY = "counts", families = list("X" = mgcv::nb(), "Y" = mgcv::nb())
    ))
})
