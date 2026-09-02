#' Test for bivariate spatial association between a point pattern and quantitative spatial omics modality in a single image
#'
#' This test is meant for the combination of a point pattern with a quantitative outcome with fixed grid location.
#' For two quantitative omics types, see \link{sbivarSingle}. For testing between point patterns, see e.g. \link[smoppix]{smoppix}
#'
#' @details If only Cx is supplied and X and Y have the same number of rows, a joint analysis is performed
#' If Cx and Ey are provided, and X and Y have the same number of rows, equality of Cx and Ey is checked.
#' If true, a joint analysis is run, with a warning.
#'
#' @param X A marked point pattern of class \link[spatstat]{ppp}
#' @param Y Matrix of omics measurements for the second, quantitative modality
#' @param Ey Coordinate matrix of dimension two, belonging to Y
#' @param method A character string, indicating which method to apply
#' @param correlation Correlation structure, passed onto \link{fitGAM}
#' @param n_points_grid,families,Gamm Passed onto \link{GAMsSingle} for the second modality fitting
#' @param wo,variogramModels,numNNs,etas,cutoff,width,returnSEsMoransI,findMaxW Parameters for the calculation of Moran's I, passed onto \link{buildWeightMat}
#' @param verbose Should info on type of analysis be printed?
#' @param normX,normY,pseudoCount Normalization parameters, passed onto \link{normMat}
#' @param featuresX,featuresY Features to be tested. Defaults to all features, but specifying them allows to test a limited feature set,
#' while using the whole matrix to calculate library sizes as offset or for normalization for Y.
#'
#' @details Ey must have rownames matching those in Y, and have two columns.
#' For GAMs, usually no normalization is needed, as the non-gaussianity is taken care of by
#' the outcome distribution, offset and link functions. Currently, identity, inverse and log-link are implemented.
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
#' @note All methods use multithreading on the cluster provided using the BiocParallel package
sbivarSinglePPP <- function(
        X, Y, Cx, Ey, method = c("Moran's I", "GAMs"),
        normX = c("none", "rel", "log"), normY = c("none", "rel", "log"), pseudoCount = 1e-8,
        etas = c(5e-6, 2e-4, 2e-2), findMaxW = FALSE, returnSEsMoransI = TRUE,
        family = gaussian(), Gamm = FALSE, featuresX = unique(marks(X, drop=FALSE)$features), featuresY = colnames(Y),
        n_points_grid = 6e2, verbose = TRUE,
        variogramModels = c("Exp", "Lin"), width = cutoff / 15, cutoff = sqrt(2) / 3,
        wo = c("Gauss", "nn"), numNNs = c(4, 8, 24),
        correlation = corGaus(form = ~ x + y, nugget = TRUE, value = c(0.9 * max(apply(Ey, 2, function(x) diff(range(x)))), 0.25))
) {
    stopifnot(
        is.numeric(n_points_grid), ncol(Cx) == 2, is.numeric(numNNs), all(numNNs > 0),
        is(family$link, "family"), family$link %in% c("identity", "log", "inverse"),
        !is.null(colnames(Y)), is.logical(Gamm), inherits(correlation, "corSpatial"),
        is.numeric(etas), all(featuresX %in% unique(marks(X, drop=FALSE)$features)),
        all(featuresY %in% colnames(Y)), !anyDuplicated(featuresX), !anyDuplicated(featuresY),
        is.logical(verbose), is.logical(findMaxW)
    )
    method <- match.arg(method)
    variogramModels <- match.arg(variogramModels, several.ok = TRUE)
    normX <- match.arg(normX)
    normY <- match.arg(normY)
    foo <- checkInputSingle(X, Y, Cx, Ey)
    Y <- normMat(Y, normY, pseudoCount)
    featuresX <- make.names(featuresX)
    featuresY <- make.names(featuresY)
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
            X = X, Y = Y, Ey = Ey, wo = wo, numNNs = selfName(numNNs),
            variogramModels = variogramModels, etas = selfName(etas), width = width,
            returnSEsMoransI = returnSEsMoransI, verbose = verbose, cutoff = cutoff, findMaxW = findMaxW,
            featuresX = featuresX, featuresY = featuresY
        ))$res
    } else if (method == "GAMs") {
        GAMsSinglePPP(
            X = X, Y = Y, Ey = Ey, families = families, n_points_grid = n_points_grid,
            verbose = verbose, featuresX = featuresX, featuresY = featuresY, Gamm = Gamm, correlation = correlation
        )
    }
    out <- cbind(out, "pAdj" = p.adjust(out[, "pVal"], method = "BH"))
    out <- addFeatureColumn(out[order(out[, "pVal"]), , drop = FALSE])
    lis <- list(
        "result" = out, "method" = method,
        "multi" = FALSE, "normX" = normX, "normY" = normY
    )
    if (method == "Moran's I") {
        lis$maxIxy <- moranRes$maxIxy
        lis$wo <- wo
        lis$wParams <- switch(wo, "Gauss" = etas, "nn" = numNNs)
    }
    if (method == "GAMs") {
        lis$families <- families
        lis$correlation <- if (Gamm) correlation
        lis$Gamm <- Gamm
    }
    return(lis)
}
