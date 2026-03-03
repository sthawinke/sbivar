#' Test for bivariate spatial association in a single image
#'
#' The tests can be applied to either disjoint or joint coordinate sets.
#'
#' @details If only Cx is supplied and X and Y have the same number of rows, a joint analysis is performed
#' If Cx and Ey are provided, and X and Y have the same number of rows, equality of Cx and Ey is checked.
#' If true, a joint analysis is run, with a warning.
#'
#' @param X,Y Matrices of omics measurements
#' @param Cx,Ey Corresponding coordinate matrices of dimension two
#' @param method A character string, indicating which method to apply
#' @param GPmethod,Quants,numLscAlts,optControl,corStruct,gpParams Passed onto \link{fitGP}
#' @param n_points_grid,families Passed onto \link{GAMsSingle}
#' @param wo,variogramModels,numNNs,etas,cutoff,width,returnSEsMoransI,findMaxW Parameters for the calculation of Moran's I, passed onto \link{buildWeightMat}
#' @param verbose Should info on type of analysis be printed?
#' @param normX,normY,pseudoCount Normalization parameters, passed onto \link{normMat}
#' @param featuresX,featuresY Features to be tested. Defaults to all features, but specifying them allows to test a limited feature set,
#' while using the whole matrix to calculate library sizes as offset or for normalization.
#'
#' @details For GAMs, usually no normalization is needed, as the non-gaussianity is taken care of by
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
#' @importFrom nlme corGaus lmeControl
#' @note All methods use multithreading on the cluster provided using the BiocParallel package
#' @seealso \link{MoransISingle}, \link{ModTtestSingle}, \link{GAMsSingle}, \link{GPsSingle}
sbivarSingle <- function(X, Y, Cx, Ey, method = c("Moran's I", "GAMs", "Modified t-test", "GPs"),
    normX = c("none", "rel", "log"), normY = c("none", "rel", "log"), pseudoCount = 1e-8,
    etas = c(5e-6, 2e-4, 2e-2), findMaxW = FALSE, returnSEsMoransI = TRUE,
    families = list("X" = gaussian(), "Y" = gaussian()), featuresX = colnames(X), featuresY = colnames(Y),
    n_points_grid = 6e2, verbose = TRUE,
    variogramModels = c("Exp", "Lin"), width = cutoff / 15, cutoff = sqrt(2) / 3,
    wo = c("Gauss", "nn"), numNNs = c(4, 8, 24),
    GPmethod = c("REML", "ML"), gpParams, Quants = c(0.005, 0.5), numLscAlts = 5,
    optControl = lmeControl(
        opt = "optim", maxIter = 5e2, msMaxIter = 5e2,
        niterEM = 1e3, msMaxEval = 1e3
    ),
    corStruct = corGaus(form = ~ x + y, nugget = TRUE, value = c(1, 0.25))) {
    stopifnot(
        is.numeric(n_points_grid), ncol(Cx) == 2, is.numeric(numNNs), all(numNNs > 0),
        all(vapply(families, FUN.VALUE = TRUE, is, "family")), is.list(optControl), !is.null(colnames(X)),
        !is.null(colnames(Y)),
        inherits(corStruct, "corGaus"), is.numeric(etas), all(featuresX %in% colnames(X)),
        all(featuresY %in% colnames(Y)), !anyDuplicated(featuresX), !anyDuplicated(featuresY),
        length(Quants) == 2, is.numeric(Quants), is.logical(verbose), is.logical(findMaxW)
    )
    if (verbose) {
        message(
            "Starting sbivar analysis of a single image on ",
            bpparam()$workers, " computing cores"
        )
    }
    foo <- checkInputSingle(X, Y, Cx, Ey)
    X <- addDimNames(X, "X")
    Y <- addDimNames(Y, "Y")
    featuresX <- make.names(featuresX)
    featuresY <- make.names(featuresY)
    method <- match.arg(method)
    variogramModels <- match.arg(variogramModels, several.ok = TRUE)
    normX <- match.arg(normX)
    normY <- match.arg(normY)
    GPmethod <- match.arg(GPmethod)
    wo <- match.arg(wo)
    if (missing(Ey)) {
        if (nrow(X) != nrow(Y)) {
            stop("Only one coordinate matrix supplied, but sample size of X and Y do not match!
                 Please supply the 'Ey' argument and consider 'Moran's I' as method.")
        }
        # Run a joint analysis
        if (verbose) {
            message("Only one coordinate matrix Cx supplied, and dimensions of X and Y do match.
                 Performing an analysis with joint coordinate sets.")
        }
        Ey <- Cx
    } else if (method == "Modified t-test") {
        stop("Two coordinate matrices supplied, whereas the modified t-test requires only one.
            If coordinates are shared, omit Ey. If coordinates are disjoint, consider using method = 'MoransI'.")
    }
    if (!identical(names(families), c("X", "Y"))) {
        stop("Name families 'X' and 'Y' for unambiguous matching")
    }
    if(is.null(rownames(Cx))){
        rownames(Cx) <- rownames(X)
    } else {
        Cx <- Cx[rownames(X),]
    }
    if(is.null(rownames(Ey))){
        rownames(Ey) <- rownames(Y)
    } else {
        Ey <- Ey[rownames(Y),]
    }
    if (!missing(gpParams)) {
        if (!is.list(gpParams) || names(gpParams) != c("X", "Y") ||
            !all(vapply(gpParams, FUN.VALUE = TRUE, function(x) {
                identical(sort(rownames(x), c("mean", "nugget", "range", "sigma")))
            })) || !all(colnames(X) %in% colnames(gpParams$X) ||
            !all(colnames(Y) %in% colnames(gpParams$Y)))) {
            stop("gpParams must be a list with names 'X' and 'Y' consisting of
            matrices with rownames 'mean', 'nugget', 'range', 'sigma',
                 and colnames the same as matrices X and Y!")
        }
    }
    if ((normX == "log" || normY == "log") && method == "GAMs") {
        warning("Normalizing data is not recommended for GAMs!
                try accounting for non-normality through the 'families' argument.", immediate. = TRUE)
    }
    colnames(Cx) <- colnames(Ey) <- c("x", "y")
    X <- normMat(X, normX, pseudoCount)
    Cx <- Cx[rownames(X), ]
    Y <- normMat(Y, normY, pseudoCount)
    Ey <- Ey[rownames(Y), ]
    out <- if (method == "Moran's I") {
        (moranRes <- MoransISingle(
            X = X[, featuresX, drop = FALSE], Y = Y[, featuresY, drop = FALSE],
            Cx = Cx, Ey = Ey, wo = wo, numNNs = selfName(numNNs),
            variogramModels = variogramModels, etas = selfName(etas), width = width,
            returnSEsMoransI = returnSEsMoransI, verbose = verbose, cutoff = cutoff, findMaxW = findMaxW,
            featuresX = featuresX, featuresY = featuresY
        ))$res
    } else if (method == "GAMs") {
        GAMsSingle(
            X = X, Y = Y, Cx = Cx, Ey = Ey, families = families, n_points_grid = n_points_grid,
            verbose = verbose, featuresX = featuresX, featuresY = featuresY
        )
    } else if (method == "GPs") {
        GPsSingle(
            X = X, Y = Y, Cx = Cx, Ey = Ey, gpParams = gpParams, Quants = Quants,
            GPmethod = GPmethod, corStruct = corStruct, optControl = optControl,
            numLscAlts = numLscAlts, verbose = verbose, featuresX = featuresX, featuresY = featuresY
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
    result <- out[order(out[, "pVal"]), ]
    lis <- list(
        "result" = result, "method" = method,
        "multi" = FALSE, "normX" = normX, "normY" = normY
    )
    if (method == "Moran's I") {
        lis$maxIxy <- moranRes$maxIxy
        lis$wo <- wo
        lis$wParams <- switch(wo,
            "Gauss" = etas,
            "nn" = numNNs
        )
    }
    if (method == "GAMs") {
        lis$families <- families
    }
    return(lis)
}
