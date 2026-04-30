#' Name a character vector after itself
#' @param x The vector to be named
#' @return The named vector
#' @export
#' @examples
#' selfName(LETTERS[1:5])
selfName <- function(x) {
    names(x) <- x
    x
}
#' Convert z-value to p-value
#'
#' @param z The z-value to be converted
#' @return The p-value
makePval <- function(z) {
    z[is.na(z)] <- 0
    tmp <- pnorm(z, lower.tail = TRUE)
    tmp[z > 0] <- pnorm(z[z > 0], lower.tail = FALSE)
    out <- 2 * unname(tmp)
    if (is.array(z)) {
        dimnames(out) <- dimnames(z)
    }
    return(out)
}
#' Scale to [0,1] range
#' @param y The vector to be scaled
#' @param na.rm passed onto \link[base]{min} and \link[base]{range}
#' @return The scaled vector
scaleZeroOne <- function(y, na.rm = TRUE) {
    (y - min(y, na.rm = na.rm)) / diff(range(y, na.rm = na.rm))
}
#' Scale to [-1,1] range
#' @inheritParams scaleZeroOne
#' @return The scaled vector
scaleMinusOne <- function(y, na.rm = TRUE) {
    scaleZeroOne(y, na.rm = na.rm) * 2 - 1
}
#' Make unique names
#' @param featX,featY vectors of feature names
#' @return A vector of names
makeNames <- function(featX, featY) {
    make.names(apply(expand.grid(featX, featY), 1, paste, collapse = "__"))
}
#' A wrapper for Matrix::bdiag maintaining names
#'
#' @param A,B Matrix to be used in \link[Matrix]{bdiag}
#' @return Same as \link[Matrix]{bdiag} but with dimnames
bdiagn <- function(A, B) {
    M <- bdiag(A, B)
    # Build new dimnames from components
    dimnames(M) <- list(c(rownames(A), rownames(B)), c(colnames(A), colnames(B)))
    M
}
#' Find trace of a matrix, of traces of an array
#'
#' A (mxm) matric has one trace (the sum of the diagonal elements), a (mxmxp) array has p traces
#'
#' @param x Matrix or array
#' @param dim Dimensions defining matrices to find traces over
#'
#' @returns A trace or vector of traces
#' @importFrom methods is
#' @importFrom Matrix diag
tr <- function(x, dim = c(1, 2)) {
    if (is.matrix(x) || is(x, "Matrix")) {
        sum(diag(x))
    } else if (is.array(x)) {
        apply(x, dim, tr)
    } else {
        stop("Trace function not implemented for ", class(x))
    }
}
#' Get all categrocial variables from a dataframe
#'
#' @param df The data frane
#'
#' @returns A character vector of variable names
getDiscreteVars <- function(df) {
    colnames(df)[!vapply(df, FUN.VALUE = TRUE, is.numeric)]
}
#' Wrapper to normalize, select feature and scale
#'
#' @details Returns vector of NA if feature not found, leading to grey in the plots
#' @param X data matrix
#' @param feat the feature name
#'
#' @returns A vector of values
scaleHelpFun <- function(X, feat) {
    if (feat %in% colnames(X)) {
        scaleZeroOne(X[, feat])
    } else {
        rep(NA, nrow(X))
    }
}
#' Normalize a data matrix, and ensure correct column names
#'
#' Normalize to relative expression, and potentially add pseudocount and log-normalize.
#' @param x The matrix
#' @param norm A character string, either "none", "log" or "rel"
#' @param pseudoCount A pseudocount added prior to log-normalization to avoid taking the log of zero
#' @details norm = "none" is pass-through, norm = "rel" divides by sample sums,
#' "log" adds a pseudocount, divides by sample sums and log-normalizes.
#' @returns A normalized matrix
#' @export
#' @examples
#' mat <- matrix(rpois(2000, lambda = 3), 40, 50)
#' nMat <- normMat(mat, norm = "rel")
normMat <- function(x, norm, pseudoCount = 1e-8) {
    stopifnot(is.matrix(x))
    if (norm == "none") {
        out <- x
    } else {
        if (any(x < 0)) {
            warning("Normalization '", norm, "'is not recommended for real valued data!")
        }
        x <- x[rowSums(x) > 0, , drop = FALSE]
        dn <- dimnames(x)
        out <- if (norm == "rel") {
            x / rowSums(x)
        } else if (norm == "log") {
            log((x + pseudoCount) / rowSums(x))
        }
        dimnames(out) <- dn
    }
    colnames(out) <- make.names(colnames(out))
    return(out)
}
#' Split a string
#'
#' @param string The string
#' @param split string to split by
#'
#' @returns A character vector of length 2
sund <- function(string, split = "__") {
    strsplit(string, split = split)[[1]]
}
#' Replace the left hand side of a formula by a fixed string
#'
#' @param x a formula
#' @param repl the replacement string
#'
#' @returns A formula
#' @importFrom stats formula
replaceLhs <- function(x, repl = "out") {
    out <- paste(repl, "~", deparse(x[[length(x)]]))
    return(formula(out))
}
#' Is there any double underscore in the character vector?
#'
#' @param charVec The character vector
#'
#' @returns A boolean
findDoubleUnderScore <- function(charVec) {
    any(grepl("__", charVec))
}
#' Extract an assay, and transpose
#'
#' @param x The SummarizedExperiment object
#' @param assayName The name of the assay
#'
#' @returns The required assay, transposed
#' @importFrom SummarizedExperiment assay
assayT <- function(x, assayName) {
    t(assay(x, assayName))
}
#' Split a SpatialExperiment object into images
#'
#' Split Spatial Experiment into a list of SpatialExperiment objects,
#' based on a variable present in the colData slot
#'
#' @param spe The SpatialExperiment object
#' @param sample_id A character vector
#'
#' @returns A list of SpatialExperiment objects
#' @importFrom SummarizedExperiment colData
splitSpatialExperiment <- function(spe, sample_id) {
    spe_list <- split(seq_len(ncol(spe)), colData(spe)[, sample_id])
    lapply(spe_list, function(cols) spe[, cols])
}
#' Print a message for the current iteration
#'
#' @param current current state of the iterator
#' @param all vector of all iterators
#'
#' @returns prints message to output
printIteration <- function(current, all) {
    message("Image ", which(all == current), " of ", length(all))
}
#' Extract data matrix
#'
#' @param X The matrix or SpatialExperiment object
#' @param assay The name of the assay
#'
#' @returns A matrix
getX <- function(X, assay) {
    if (inherits(X, "SpatialExperiment")) {
        assayT(X, assay)
    } else if (is.list(X)) {
        lapply(X, getX, assay = assay)
    } else {
        X
    }
}
#' Extract coordinate matrix
#'
#' @param X The matrix or SpatialExperiment object
#' @param Cx The coordinate matrix
#'
#' @returns A coordinate matrix
getSpatialCoords <- function(X, Cx) {
    out <- if (inherits(X, "SpatialExperiment")) {
        SpatialExperiment::spatialCoords(X)
    } else if (is.list(X)) {
        lapply(selfName(names(X)), function(i) {
            getSpatialCoords(X[[i]], Cx[[i]])
        })
    } else {
        Cx
    }
    if (is.matrix(out)) {
        rownames(out) <- rownames(X)
    }
    out
}
#' Move two sets of coordinates to 0-1. without shifting them with respect to each other
#'
#' @param Cx,Ey Coordinate matrices
#'
#' @returns A list with the shifted and scaled coordinates
moveTwoCoords <- function(Cx, Ey) {
    # Bring coordinates to 0-1
    minX <- min(c(Cx[, "x"], Ey[, "x"]))
    minY <- min(c(Cx[, "y"], Ey[, "y"]))
    Cx[, "x"] <- Cx[, "x"] - minX
    Cx[, "y"] <- Cx[, "y"] - minY
    Ey[, "x"] <- Ey[, "x"] - minX
    Ey[, "y"] <- Ey[, "y"] - minY
    MaxCoord <- max(c(Cx, Ey))
    Cx <- Cx / MaxCoord
    Ey <- Ey / MaxCoord
    list("Cx" = Cx, "Ey" = Ey)
}
#' Print feature progress
#'
#' @param feat The current feature
#' @param allFeats Vector of all features
#' @param verbose Boolean, should output be printed?

#' @returns Prints progress message to output
printProgress <- function(feat, allFeats, verbose) {
    if (verbose && (la <- length(allFeats)) >= 10) {
        keyFeats <- allFeats[round(la * seq(0.1, 1, by = 0.1))]
        if (!is.na(id <- match(feat, keyFeats))) {
            message(round(10 * id), "% of tests completed")
        }
    }
}
#' Add columns with feature names to a matrix by splitting rownames, and remove rownames
#'
#' @param x A results matrix
#'
#' @returns The matrix extended with two columns of feature names in front
addFeatureColumn <- function(x) {
    featureMet <- t(vapply(rownames(x), FUN.VALUE = character(2), sund))
    dimnames(featureMet) <- list(NULL, c("Modality_X", "Modality_Y"))
    rownames(x) <- NULL
    data.frame(featureMet, x)
}
