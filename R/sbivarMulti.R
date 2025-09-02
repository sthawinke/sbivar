#' Apply one or more tests for bivariate spatial association
#'
#' @param Xl,Yl Lists of matrices of omics measurements
#' @param Cxl,Eyl Lists of corresponding coordinate matrices of dimension two
#' @param method A character string, indicating which method to apply
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
#' ims = 8
#' X = lapply(selfName(ims), function(i){n = rpois(1, n)
#'  matrix(rnorm(n*p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
#' })
#' Y = lapply(selfName(ims), function(i){m = rpois(1, m)
#'  matrix(rnorm(m*k), m, k, dimnames = list(NULL, paste0("Y", seq_len(k))))
#' })
#' Cx = lapply(X, function(x){n = nrow(x)
#'     matrix(runif(n*2), n, 2, dimnames = list(NULL, c("x", "y")))
#' })
#' Ey = lapply(Y, function(y){m = nrow(y)
#'     matrix(runif(m*2), m, 2, dimnames = list(NULL, c("x", "y")))
#' })
#' resGAMs = sbivarMulti(X, Y, Cx, Ey, method = "GAMs")
#' resMoran = sbivarMulti(X, Y, Cx, Ey, method = "Moran")
sbivarMulti = function(Xl, Yl, Cxl, Eyl, method = c("GAMs", "Correlation", "GPs", "Moran's I"),
                        wo = c("distance", "nn"), numNN = 8, n_points_grid = 5e2,
                        families = list("X" = gaussian(), "Y" = gaussian())){
    method = match.arg(method)
    wo = match.arg(wo)
    stopifnot(is.numeric(numNN), is.numeric(n_points_grid), is.character(wo),
              is.character(method), all(vapply(families, FUN.VALUE = TRUE, is, "family")))
    if(length(Xl)!=length(Cxl)){
        stop("Length of outcome matrices Xl and their coordinates Cxl do not match!")
    }
    if(length(Yl)!=length(Eyl)){
        stop("Length of outcome matrices Yl and their coordinates Eyl do not match!")
    }
    if(!all(vapply(Xl, nrow, FUN.VALUE = 0) == vapply(Cxl, nrow, FUN.VALUE = 0))){
        stop("Sample size of matrices Xl and their coordinates Cxl do not match!")
    }
    if(!all(vapply(Yl, nrow, FUN.VALUE = 0) == vapply(Eyl, nrow, FUN.VALUE = 0))){
        stop("Sample size of matrices Yl and their coordinates Eyl do not match!")
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
    out = if(method == "GAMs"){
        # wrapGAMs(X = X, Y = Y, Cx = Cx, Ey = Ey, families = families,
        #          n_points_grid = n_points_grid)
    } else if(method == "GPs"){

    } else if(method == "Modified t-test"){
 #       wrapModTtest(X = X, Y = Y, Cx = Cx, Ey = Ey, mapToFinest = FALSE)
    } else if (method =="Moran's I"){
    }
    #out[order(out[, "pVal"]),]
}

