#' Estimate measures of bivariate spatial association for multiple images
#'
#' This function calculates measures of spatial association for every image.
#' The resulting estimates can then be analysed further using the \link{fitLinModels} function.
#'
#' @param Xl,Yl Lists of matrices of omics measurements
#' @param Cxl,Eyl Lists of corresponding coordinate matrices of dimension two
#' @param method A character string, indicating which method to apply
#' @param families,n_points_grid Passed onto \link{GAMsMulti}
#' @param wo,numNNs,etas Passed onto \link{MoransISingle}
#' @returns A list containing
#' \item{estimates}{The estimated measures of association}
#' \item{multi}{TRUE, a flag for the type of analysis}
#' \item{normX,normY,method}{As provided}
#' \item{families,wo,wParams}{Optional, as provided. wParams are either etas or numNNs}
#' @note All methods use multithreading on the cluster provided using the BiocParallel package
#' @inheritParams sbivarSingle
#' @inheritParams buildWeightMat
#' @inheritParams MoransISingle
#' @importFrom BiocParallel bpparam
#' @seealso \link{fitLinModels}, \link{MoransIMulti}, \link{correlationsMulti}, \link{GAMsMulti}
sbivarMulti <- function(
      Xl, Yl, Cxl, Eyl, families = list("X" = gaussian(), "Y" = gaussian()),
      method = c("Moran's I", "GAMs", "Correlation"), wo = c("Gauss", "nn"),
      numNNs = c(4, 8, 24), etas = c(5e-6, 2e-4, 2e-2),
      normX = c("none", "rel", "log"), normY = c("none", "rel", "log"),
      variogramModels = c("Exp", "Lin"), width = cutoff / 15, cutoff = sqrt(2) / 3,
      pseudoCount = 1e-8, n_points_grid = 6e2, verbose = TRUE, findVariances = FALSE, findMaxW = TRUE
) {
    method <- match.arg(method)
    wo <- match.arg(wo)
    normX <- match.arg(normX)
    normY <- match.arg(normY)
    stopifnot(
        is.numeric(numNNs), all(numNNs > 0), is.numeric(etas), is.numeric(n_points_grid), is.logical(verbose),
        is.character(method), all(vapply(families, FUN.VALUE = TRUE, is, "family"))
    )
    Xl <- lapply(Xl, addDimNames, "X")
    Yl <- lapply(Yl, addDimNames, "Y")
    foo <- checkInputMulti(Xl, Yl, Cxl, Eyl, checkCoords = ccs <- (method != "Correlation"))
    if (ccs) {
        if (missing(Eyl)) {
            message("Only one coordinate matrix list Cxl supplied, and dimensions of X and Y do match.
                 Performing an analysis with joint coordinate sets.")
            Eyl <- Cxl
        }
        Cxl <- mapply(Cxl, Xl, SIMPLIFY = FALSE, FUN = tmpFun <- function(cx, x) {
            colnames(cx) <- c("x", "y")
            rownames(cx) <- rownames(x)
            cx
        })
        Eyl <- mapply(Eyl, Yl, SIMPLIFY = FALSE, FUN = tmpFun)
    } else if (!missing(Cxl)) {
        warning("Correlation analysis will ignore coordinate matrices provided.
                Consider providing another 'method' argument for a full spatial analysis", immediate. = TRUE)
    }
    if (verbose) {
        message(
            "Starting sbivar analysis (", method, ") of ", length(Xl), " images on ",
            bpparam()$workers, " computing cores"
        )
    }
    # Normalization
    Xl <- lapply(Xl, normMat, normX, pseudoCount)
    Yl <- lapply(Yl, normMat, normY, pseudoCount)
    if (ccs) {
        Cxl <- lapply(selfName(names(Cxl)), function(i) {
            Cxl[[i]][rownames(Xl[[i]]), ]
        })
        Eyl <- lapply(selfName(names(Eyl)), function(i) {
            Eyl[[i]][rownames(Yl[[i]]), ]
        })
    }
    out <- if (method == "Moran's I") {
        MoransIMulti(Xl, Yl, Cxl, Eyl,
            wo = wo, numNNs = numNNs, verbose = verbose, findVariances = findVariances, findMaxW = findMaxW,
            etas = etas, variogramModels = variogramModels, width = width, cutoff = cutoff
        )
    } else if (method == "GAMs") {
        GAMsMulti(Xl, Yl, Cxl, Eyl,
            families = families, findVariances = findVariances,
            n_points_grid = n_points_grid, verbose = verbose
        )
    } else if (method == "Correlation") {
        correlationsMulti(Xl, Yl, verbose = verbose)
    }
    out <- list(
        "estimates" = out, "method" = method, "multi" = TRUE,
        "normX" = normX, "normY" = normY
    )
    if (method == "GAMs") {
        out$families <- families
    } else if (method == "Moran's I") {
        out$wo <- wo
        out$wParams <- selfName(switch(wo,
            "Gauss" = etas,
            "nn" = numNNs
        ))
        out$returnSEsMoransI <- findVariances
    }
    return(out)
}
