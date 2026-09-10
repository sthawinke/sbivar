#' Test for bivariate spatial association between a point pattern and quantitative spatial omics modality in a single image
#'
#' This test is meant for the combination of a point pattern with a quantitative outcome with fixed grid location.
#' For two quantitative omics types, see \link{sbivarSingle}. For testing between point patterns, see e.g. \link[smoppix]{smoppix}
#'
#' @param X A marked point pattern of class \link[spatstat.geom]{ppp}
#' @param Y Matrix of omics measurements for the second, quantitative modality
#' @param Ey Coordinate matrix of dimension two, belonging to Y
#' @param method A character string, indicating which method to apply
#' @param correlation Correlation structure, passed onto \link{fitGAM}
#' @param n_points_grid,families,Gamm Passed onto \link{GAMsSingle} for the second modality fitting
#' @param wo,variogramModels,etas,cutoff,width,returnSEsMoransI,findVariances Parameters for the calculation of Moran's I, passed onto \link{buildWeightMat}
#' @param verbose Should info on type of analysis be printed?
#' @param normY,pseudoCount Normalization parameters, passed onto \link{normMat}
#' @param featuresX,featuresY Features to be tested. Defaults to all features, but specifying them allows to test a limited feature set,
#' while using the whole matrix to calculate library sizes as offset or for normalization for Y.
#'
#' @details Ey must have rownames matching those in Y, and have two columns.
#' For GAMs, usually no normalization is needed, as the non-gaussianity is taken care of by
#' the outcome distribution, offset and link functions. Currently, identity, inverse and log-link are implemented.
#' For constructing weight matrices, only the Gaussian weights are currently implemented, as nearest-neighbour weights are not suitable for point patterns.
#' For this purpose, eta values larger than 2e-3 are not meaningful as they lead to weight matrices covering the whole measurement area, to whicht the single molecules are restricted by design.
#' As this may lead to false findings, a warning is thrown.
#'
#' @returns A list with at least the following components
#' \item{result}{A matrix which contains at least a p-values ("pVal") and a
#' Benjamini-Hochberg adjusted p-value ("pAdj"), sorted by increasing p-value.}
#' \item{multi}{FALSE, a flag for the type of analysis}
#' \item{method,normX,normY}{As provided}
#' \item{families,wo,wParams}{Optional, as provided. wParams are either etas or numNNs}
#' @importFrom stats p.adjust
#' @importFrom methods is
#' @importFrom nlme corRatio corGaus corSpher corExp corLin lmeControl
#' @importFrom BiocParallel bpparam bpworkers
#' @importFrom spatstat.geom marks<-
#' @note All methods use multithreading on the cluster provided using the BiocParallel package
sbivarSinglePPP <- function(
      X, Y, Ey, method = c("Moran's I", "GAMs"),
      normY = c("none", "rel", "log"), pseudoCount = 1e-8,
      etas = c(5e-6, 4e-5, 2e-4), returnSEsMoransI = TRUE, findVariances = TRUE,
      families = list("Y" = gaussian()), Gamm = FALSE, featuresX = getFeaturesX(X), featuresY = colnames(Y),
      n_points_grid = 6e2, verbose = TRUE, wo = "Gauss",
      variogramModels = c("Exp", "Lin"), width = cutoff / 15, cutoff = sqrt(2) / 3,
      correlation = corGaus(form = ~ x + y, nugget = TRUE, value = c(0.9 * max(apply(Ey, 2, function(x) diff(range(x)))), 0.25))
) {
    stopifnot(
        is.numeric(n_points_grid),
        names(families) == "Y",
        is(families[["Y"]], "family"), families[["Y"]]$link %in% c("identity", "log", "inverse"),
        !is.null(colnames(Y)), is.logical(Gamm), inherits(correlation, "corSpatial"),
        is.numeric(etas), all(featuresX %in% getFeaturesX(X)),
        all(featuresY %in% colnames(Y)), !anyDuplicated(featuresX), !anyDuplicated(featuresY),
        is.logical(verbose)
    )
    if (any(etas > 2e-3)) {
        warning("Eta values larger than 2e-3 are not meaningful for point patterns and may lead to false positive findings!")
    }
    method <- match.arg(method)
    variogramModels <- match.arg(variogramModels, several.ok = TRUE)
    normY <- match.arg(normY)
    foo <- checkInputSingle(X, Y, Cx = NULL, Ey)
    Y <- normMat(Y, normY, pseudoCount)
    featuresX <- make.names(featuresX)
    marks(X, drop = FALSE)$feature <- make.names(marks(X, drop = FALSE)$feature)
    featuresY <- make.names(featuresY)
    Ey <- Ey[rownames(Y), ]
    wo <- match.arg(wo)
    if (((normY %in% c("rel", "log"))) && method == "GAMs") {
        warning("Normalizing data is not recommended for GAMs!
                try accounting for non-normality through the 'families' argument.", immediate. = TRUE)
    }
    if (verbose) {
        message(
            "Starting sbivar analysis on point pattern + quantitative outcome of a single image on ",
            bpworkers(bpparam()), " computing cores"
        )
    }
    out <- if (method == "Moran's I") {
        (moranRes <- MoransISinglePPP(
            X = X, Y = Y, Ey = Ey, wo = wo,
            variogramModels = variogramModels, etas = selfName(etas), width = width,
            returnSEsMoransI = returnSEsMoransI, verbose = verbose, cutoff = cutoff,
            featuresX = featuresX, featuresY = featuresY, findVariances = findVariances
        ))$res
    } else if (method == "GAMs") {
        GAMsSingle(
            X = X, Y = Y, Ey = Ey, families = families, n_points_grid = n_points_grid,
            verbose = verbose, featuresX = featuresX, featuresY = featuresY, Gamm = Gamm, correlation = correlation
        )
    }
    out <- cbind(out, "pAdj" = p.adjust(out[, "pVal"], method = "BH"))
    out <- addFeatureColumn(out[order(out[, "pVal"]), , drop = FALSE])
    lis <- list(
        "result" = out, "method" = method,
        "multi" = FALSE, "normY" = normY
    )
    if (method == "Moran's I") {
        lis$wo <- wo
        lis$wParams <- switch(wo,
            "Gauss" = etas
        )
    }
    if (method == "GAMs") {
        lis$families <- families
        lis$correlation <- if (Gamm) correlation
        lis$Gamm <- Gamm
    }
    return(lis)
}
