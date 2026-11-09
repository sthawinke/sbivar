#' Test for bivariate spatial association in a single image
#'
#' The tests can be applied to either joint or disjoint coordinate sets.
#'
#' @details If only Cx is supplied and X and Y have the same number of rows, a joint analysis is performed
#' If Cx and Ey are provided, and X and Y have the same number of rows, equality of Cx and Ey is checked.
#' If true, a joint analysis is run, with a warning.
#'
#' @param X,Y Matrices of omics measurements
#' @param Cx,Ey Coordinate matrices of dimension two, belonging to X and Y respectively
#' @param method A character string, indicating which method to apply
#' @param correlation Correlation structure, passed onto \link{fitGP} or \link{fitGAM}
#' @param verbose Should info on type of analysis be printed?
#' @param normX,normY,pseudoCount Normalization parameters, passed onto \link{normMat}
#' @param featuresX,featuresY Features to be tested. Defaults to all features, but specifying them allows to test a limited feature set,
#' while using the whole matrix to calculate library sizes as offset or for normalization.
#'
#' @details X and Y need to have rownames for matching to the coordinates, and column names for identifying the features.
#' Cx and Ey must have rownames matching those in X and Y, and have two columns.
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
#' @importFrom methods is new
#' @importFrom nlme corRatio corGaus corSpher corExp corLin lmeControl
#' @importFrom BiocParallel bpparam bpworkers
#' @note All methods use multithreading on the cluster provided using the BiocParallel package.
#' @seealso \link{MoransISingle}, \link{ModTtestSingle}, \link{GAMsSingle}, \link{GPsSingle}
sbivarSingle <- function(
      X, Y, Cx, Ey, method = c("Moran's I", "GAMs", "Modified t-test", "GPs"),
      normX = c("none", "rel", "log"), normY = c("none", "rel", "log"), pseudoCount = 1e-8,
      moranList = moranControl(), GPlist = GPcontrol(Cx = Cx), GAMlist = GAMcontrol(Cx = Cx), gpParams,
      featuresX = colnames(X), featuresY = colnames(Y),
      verbose = TRUE
) {
    method <- match.arg(method)
    variogramModels <- match.arg(moranList$variogramModels, choices = c("Exp", "Lin"), several.ok = TRUE)
    normX <- match.arg(normX)
    normY <- match.arg(normY)
    stopifnot(
        ncol(Cx) == 2, !is.null(colnames(Y)), all(featuresX %in% colnames(X)),
        all(featuresY %in% colnames(Y)), !anyDuplicated(featuresX), !anyDuplicated(featuresY),
        !is.null(colnames(X)), is.logical(verbose),
        is.numeric(moranList$numNNs), all(moranList$numNNs > 0), is.numeric(moranList$etas),
        is.logical(moranList$findMaxW),
        is.numeric(GAMlist$n_points_grid),
        all(vapply(GAMlist$families, FUN.VALUE = character(1), function(x) x$link) %in% c("identity", "log", "inverse")),
        all(vapply(GAMlist$families, FUN.VALUE = TRUE, is, "family")),
        is.logical(GAMlist$Gamm), inherits(GAMlist$correlation, "corSpatial"),
        inherits(GPlist$correlation, "corGaus"), is.list(GPlist$optControl),
        length(GPlist$Quants) == 2, is.numeric(GPlist$Quants)
    )
    invisible(checkInputSingle(X, Y, Cx, Ey))
    colnames(Cx) <- c("x", "y")
    GPmethod <- match.arg(GPlist$GPmethod, choices = c("ML", "REML"))
    X <- normMat(X, normX, pseudoCount)
    Y <- normMat(Y, normY, pseudoCount)
    Cx <- Cx[rownames(X), ] # Match coordinates to observations through rownames
    featuresX <- make.names(featuresX)
    featuresY <- make.names(featuresY)
    wo <- match.arg(wo)
    if (missing(Ey)) {
        if (nrow(X) != nrow(Y)) {
            stop(strwrapPaste("Only one coordinate matrix supplied, but sample size of X and Y do not match!
                 Please supply the 'Ey' argument and consider 'Moran's I' as method."))
        }
        # Run a joint analysis
        if (verbose) {
            message(strwrapPaste("Only one coordinate matrix Cx supplied, and dimensions of X and Y do match.
                 Performing an analysis with joint coordinate sets."))
        }
        Ey <- Cx
    } else if (method == "Modified t-test") {
        stop(strwrapPaste("Two coordinate matrices supplied, whereas the modified t-test requires only one.
            If coordinates are shared, omit Ey. If coordinates are disjoint, consider using method = 'MoransI'."))
    } else {
        Ey <- Ey[rownames(Y), ]
        colnames(Ey) <- c("x", "y")
    }
    if (!identical(names(families), c("X", "Y"))) {
        stop("Name families 'X' and 'Y' for unambiguous matching")
    }
    if (!missing(gpParams)) {
        if (!is.list(gpParams) || !identical(names(gpParams), c("X", "Y")) ||
            !all(vapply(gpParams, FUN.VALUE = TRUE, function(x) {
                setequal(rownames(x), c("mean", "nugget", "range", "sigma"))
            })) || !all(colnames(X) %in% colnames(gpParams$X) ||
            !all(colnames(Y) %in% colnames(gpParams$Y)))) {
            stop(strwrapPaste("gpParams must be a list with names 'X' and 'Y' consisting of
            matrices with rownames 'mean', 'nugget', 'range', 'sigma',
                 and colnames the same as matrices X and Y!"))
        }
    }
    if (((normX %in% c("rel", "log")) || (normY %in% c("rel", "log"))) && method == "GAMs") {
        warning(strwrapPaste("Normalizing data is not recommended for GAMs!
                try accounting for non-normality through the 'families' argument."), immediate. = TRUE)
    }
    if (verbose) {
        message(
            strwrapPaste(
                "Starting sbivar analysis of a single image on ",
                bpworkers(bpparam()), " computing cores"
            )
        )
    }
    out <- if (method == "Moran's I") {
        (moranRes <- MoransISingle(
            X = X, Y = Y, Cx = Cx, Ey = Ey, wo = moranList$wo, numNNs = selfName(moranList$numNNs),
            variogramModels = moranList$variogramModels, etas = selfName(moranList$etas), width = moranList$width,
            returnSEsMoransI = moranList$returnSEsMoransI, verbose = verbose, cutoff = moranList$cutoff, findMaxW = moranList$findMaxW,
            featuresX = featuresX, featuresY = featuresY
        ))$res
    } else if (method == "GAMs") {
        GAMsSingle(
            X = X, Y = Y, Cx = Cx, Ey = Ey, families = GAMlist$families, n_points_grid = GAMlist$n_points_grid,
            verbose = verbose, featuresX = featuresX, featuresY = featuresY, Gamm = GAMlist$Gamm, correlation = GAMlist$correlation
        )
    } else if (method == "GPs") {
        GPsSingle(
            X = X, Y = Y, Cx = Cx, Ey = Ey, gpParams = gpParams, Quants = GPlist$Quants,
            GPmethod = GPlist$GPmethod, correlation = GPlist$correlation, optControl = GPlist$optControl,
            numLscAlts = GPlist$numLscAlts, verbose = verbose, featuresX = featuresX, featuresY = featuresY
        )
    } else if (method == "Modified t-test") {
        sharedNames <- intersect(rownames(X), rownames(Y))
        ModTtestSingle(
            X = X[sharedNames, featuresX, drop = FALSE],
            Y = Y[sharedNames, featuresY, drop = FALSE],
            Cx = Cx[sharedNames, ], verbose = verbose
        )
    }
    out <- cbind(out, "pAdj" = p.adjust(out[, "pVal"], method = "BH"))
    out <- addFeatureColumn(out[order(out[, "pVal"]), , drop = FALSE])
    res <- if (method == "Moran's I") {
        new("sbivarResultsMoransI",
            "result" = out, "method" = method, "normX" = normX, "multi" = FALSE,
            "normY" = normY, "maxIxy" = moranRes$maxIxy, "wo" = wo, "returnSEsMoransI" = moranList$returnSEsMoransI,
            "wParams" = switch(wo,
                "Gauss" = moranList$etas,
                "nn" = moranList$numNNs
            )
        )
    } else if (method == "GAMs") {
        new("sbivarResultsGAMs",
            "result" = out, "method" = method, "multi" = FALSE,
            "normX" = normX, "normY" = normY, "families" = GAMlist$families,
            "correlation" = if (Gamm) GAMlist$correlation, "Gamm" = GAMlist$Gamm
        )
    } else {
        new("sbivarResults",
            "result" = out, "method" = method, "multi" = FALSE,
            "normX" = normX, "normY" = normY
        )
    }
    return(res)
}
