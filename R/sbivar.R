#' Apply one or more tests for bivariate spatial association
#'
#' @param X,Y Matrices of omics measurements
#' @param Cx,Ey Corresponding coordinate matrices of dimension two
#' @param method A character string, indicating which method to apply
#' @param wMat Optional, a weight matrix for calculating Moran's I
#' @param wo,numNN Passed on to \link{buildWeightMat}.
#' @param families A vector of length 2 giving outcome values.
#' @param n_points_grid The number of points in the new grid for the GAMs to be
#' evaluated on.
#'
#' @returns A list containing at least an entry "result", which contains p-values,
#' adjusted p-values and measures of association, sorted by increasing p-value.
#' @export
#'
#' @examples
#' n=1e2;m=2e2;p=10;k=12
#' X = matrix(rnorm(n*p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
#' Y = matrix(rnorm(m*k), m, k, dimnames = list(NULL, paste0("Y", seq_len(k))))
#' Cx = matrix(runif(n*2), n, 2)
#' Ey = matrix(runif(m*2), m, 2)
#' colnames(Cx) = colnames(Ey) = c("x", "y")
#' resGAMs = sbivar(X, Y, Cx, Ey, method = "GAMs")
sbivar = function(X, Y, Cx, Ey, method = c("GAMs", "Modified t-test", "GPs"),
                  wMat, wo = "distance", numNN = 8, n_points_grid = 5e2,
                  families = list("X" = gaussian(), "Y" = gaussian())){
    stopifnot(is.numeric(numNN), is.numeric(n_points_grid), is.character(wo),
              is.character(method), all(vapply(families, FUN.VALUE = TRUE, is, "family")))
    n = nrow(X);m = nrow(Y);p = ncol(X);k=ncol(Y)
    if(n!=nrow(Cx)){
        stop("Dimensions of X and its coordinates Cx do not match!")
    }
    if(m!=nrow(Ey)){
        stop("Dimensions of Y and its coordinates Ey do not match!")
    }
    if(ncol(Cx)!=2 || ncol(Ey)!=2){
        stop("Coordinate matrices must be of dimension 2!")
    }
    if(!identical(names(families), c("X", "Y"))){
        stop("Name families 'X' and 'Y' for unambiguous matching")
    }
    if(is.null(colnames(X))){
        colnames(X) = paste0("X", seq_len(p))
    }
    if(is.null(colnames(Y))){
        colnames(Y) = paste0("Y", seq_len(k))
    }
    colnames(Cx) = colnames(Ey) = c("x", "y")
    method = match.arg(method)
    out = if(method == "Moran's I"){

    } else if(method == "GAMs"){
        wrapGAMs(X = X, Y = Y, Cx = Cx, Ey = Ey, families = families,
                 n_points_grid = n_points_grid)
    } else if(method == "GPs"){

    } else if(method == "Modified t-test"){
        wrapModTtest(X = X, Y = Y, Cx = Cx, Ey = Ey, mapToFinest = FALSE)
    }
}

