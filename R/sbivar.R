#' Apply one or more tests for bivariate spatial association
#'
#' @param X,Y Matrices of omics measurements
#' @param Cx,Ey Corresponding coordinate matrices of dimension two
#' @param method A character string, indicating which method to apply
#' @param wMat Optional, a weight matrix for calculating Moran's I
#' @param wo,numNN Passed on to \link{buildWeightMat}.
#' @param n_points_grid
#'
#' @returns A list containing at least an entry "result", which contains p-values,
#' adjusted p-values and measures of association, sorted by increasing p-value.
#' @export
#'
#' @examples
sbivar = function(X, Y, Cx, Ey, method = c("GAMs", "Modified t-test", "GPs"), wMat, wo = "distance", numNN = 8,
                  n_points_grid = 5e2){
    stopifnot(is.numeric(numNN), is.numeric(n_points_grid), is.character(wo),
              is.character(methods))
    n = nrow(X);m = nrow(Y);p = ncol(X);k=ncol(y)
    if(n!=nrow(Cx)){
        stop("Dimensions of X and its coordinates Cx do not match!")
    }
    if(m!=nrow(Ey)){
        stop("Dimensions of Y and its coordinates Ey do not match!")
    }
    if(ncol(Cx)!=2 || ncol(Ey)!=2){
        stop("Coordinate matrices must be of dimension 2!")
    }
    method = match.arg(methods)
    out = if(method == "Moran's I"){

    } else if(method == "GAMs"){

    } else if(method == "GPs"){

    } else if(method == "Modified t-test"){

    }
}

