#' Check input of matrices and lists
#'
#' Checks dimensions of matrices and lists, to be used for all exported functions
#'
#' @inheritParams sbivarSingle
#' @inheritParams sbivarMulti
#'
#' @returns Throws error when a fault is found with the input, otherwise finishes silently
checkInputSingle <- function(X, Y, Cx, Ey) {
    n <- nrow(X)
    m <- nrow(Y)
    if (n != nrow(Cx)) {
        stop("Dimensions of X and its coordinates Cx do not match!")
    }
    if (ncol(Cx) != 2) {
        stop("Coordinate matrices must be of dimension 2!")
    }
    if (is.null(colnames(X))) {
        stop("Feature matrix X lacks column names!")
    }
    if (is.null(colnames(Y))) {
        stop("Feature matrix Y lacks column names!")
    }
    if (is.null(rownames(X))) {
        stop("Feature matrix X lacks column names!")
    }
    if (is.null(rownames(Y))) {
        stop("Feature matrix Y lacks column names!")
    }
    if (missing(Ey)) {
        if (n != m) {
            stop("Only one coordinate matrix Cx supplied, and dimensions of X and Y do not match.
                 Please provide the coordinates of Y too through the Ey argument.")
        }
    } else {
        if (m != nrow(Ey)) {
            stop("Dimensions of Y and its coordinates Ey do not match!")
        }
        if (ncol(Ey) != 2) {
            stop("Coordinate matrices must be of dimension 2!")
        }
        if (!identical(sort(rownames(Y)), sort(rownames(Ey)))) {
            stop("Rownames of Y and Ey do not match")
        }
    }
    if (findDoubleUnderScore(c(colnames(X), colnames(Y)))) {
        stop("Double underscores found in feature names. Please change the names,
             as the double underscore is used in this package to separate feature pairs!")
    }
    if (!identical(sort(rownames(X)), sort(rownames(Cx)))) {
        stop("Rownames of X and Cx do not match")
    }
}
checkInputMulti <- function(Xl, Yl, Cxl, Eyl, checkCoords = TRUE) {
    if (length(Xl) == 1) {
        stop("Lists of length 1 not allowed, please convert to matrix!")
    }
    if (any(vapply(Xl, FUN.VALUE = FALSE, function(x) is.null(colnames(x))))) {
        stop("Some feature names are missing in Xl!")
    }
    if (any(vapply(Yl, FUN.VALUE = FALSE, function(x) is.null(colnames(x))))) {
        stop("Some feature names are missing in Yl!")
    }
    if (any(vapply(Xl, FUN.VALUE = FALSE, function(x) is.null(rownames(x))))) {
        stop("Some sample names are missing in Xl!")
    }
    if (any(vapply(Yl, FUN.VALUE = FALSE, function(x) is.null(rownames(x))))) {
        stop("Some sample names are missing in Yl!")
    }
    if (findDoubleUnderScore(unlist(lapply(c(Xl, Yl), colnames)))) {
        stop("Double underscores found in feature names. Please change the names,
             as the double underscore is used in this package to separate feature pairs!")
    }
    if (!all(vapply(Xl, FUN.VALUE = 0L, nrow) == vapply(Yl, FUN.VALUE = 0L, nrow))) {
        if (missing(Eyl)) {
            stop("If Eyl is not supplied, all elements of Xl and Yl must have the number of samples!")
        } else if (!checkCoords) {
            stop("For Pearson correlation analysis, all elements of Xl and Yl must have the number of samples!")
        }
    }
    if (checkCoords) {
        if (length(Xl) != length(Cxl)) {
            stop("Length of outcome matrices Xl and their coordinates Cxl do not match!")
        }
        if (!all(vapply(Xl, nrow, FUN.VALUE = 0) == vapply(Cxl, nrow, FUN.VALUE = 0))) {
            stop("Sample size of matrices Xl and their coordinates Cxl do not match!")
        }
        if (is.null(names(Xl)) || is.null(names(Yl)) || is.null(names(Cxl))) {
            stop("Xl, Yl and Cxl must be named lists")
        }
        if (!identical(names(Xl), names(Yl)) || !identical(names(Yl), names(Cxl))) {
            stop("All names of Xl, Yl and Cxl must be identical")
        }
        if (any(!mapply(Xl, Cxl, FUN = function(x, y) {
            identical(sort(rownames(x)), sort(rownames(y)))
        }))) {
            stop("Not all sample names are identical in Xl and Cxl!")
        }
        if (!missing(Eyl)) {
            if (length(Yl) != length(Eyl)) {
                stop("Length of outcome matrices Yl and their coordinates Eyl do not match!")
            }
            if (!all(vapply(Yl, nrow, FUN.VALUE = 0) == vapply(Eyl, nrow, FUN.VALUE = 0))) {
                stop("Sample size of matrices Yl and their coordinates Eyl do not match!")
            }
            if (is.null(names(Eyl))) {
                stop("Eyl must be a named list")
            }
            if (!missing(Eyl) && !identical(names(Eyl), names(Xl))) {
                stop("Eyl must be named identically to Xl, Yl and Cxl")
            }
            if (any(!mapply(Yl, Eyl, FUN = function(x, y) {
                identical(sort(rownames(x)), sort(rownames(y)))
            }))) {
                stop("Not all sample names are identical in Yl and Eyl!")
            }
        }
    }
}
