#' Check input of matrices and lists
#'
#' Checks dimensions of matrices and lists, to be used for all exported functions
#'
#' @inheritParams sbivarSingle
#' @inheritParams sbivarMulti
#'
#' @returns Throws error when a fault is found with the input, otherwise finishes silently
checkInputSingle = function(X, Y, Cx, Ey){
    n = nrow(X);m = nrow(Y)
    if(n!=nrow(Cx)){
        stop("Dimensions of X and its coordinates Cx do not match!")
    }
    if(m!=nrow(Ey)){
        stop("Dimensions of Y and its coordinates Ey do not match!")
    }
    if(ncol(Cx)!=2 || ncol(Ey)!=2){
        stop("Coordinate matrices must be of dimension 2!")
    }
    if(n!=m){
        stop("Only one coordinate matrix Cx supplied, and dimensions of X and Y do not match.
                 Please provide the coordinates of Y too through the Ey argument.")
    }
}
checkInputMulti = function(Xl, Yl, Cxl, Eyl){
    if(length(Xl)!=length(Cxl)){
        stop("Length of outcome matrices Xl and their coordinates Cxl do not match!")
    }
    if(!all(vapply(Xl, nrow, FUN.VALUE = 0) == vapply(Cxl, nrow, FUN.VALUE = 0))){
        stop("Sample size of matrices Xl and their coordinates Cxl do not match!")
    }
    if(missing(Eyl)){
        if(!all(vapply(Xl, FUN.VALUE = 0L, nrow) == vapply(Yl, FUN.VALUE = 0L, nrow))){
            stop("If Eyl is not supplied, all elements of Xl and Yl must have the number of samples!")
        }
    } else {
        if(length(Yl)!=length(Eyl)){
            stop("Length of outcome matrices Yl and their coordinates Eyl do not match!")
        }
        if(!all(vapply(Yl, nrow, FUN.VALUE = 0) == vapply(Eyl, nrow, FUN.VALUE = 0))){
            stop("Sample size of matrices Yl and their coordinates Eyl do not match!")
        }
    }
}
