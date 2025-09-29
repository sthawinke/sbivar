#' Apply one or more tests for bivariate spatial association
#'
#' This function calculates measures of spatial association for every image.
#' The resulting estimates can then be analysed further usign the \link{fitLinModels} function.
#'
#' @param Xl,Yl Lists of matrices of omics measurements
#' @param Cxl,Eyl Lists of corresponding coordinate matrices of dimension two
#' @param method A character string, indicating which method to apply
#' @param families A vector of length 2 giving outcome values.
#' @param n_points_grid The number of points in the new grid for the GAMs to be
#' evaluated on.
#'
#' @returns A list containing
#' \item{estimates}{The estimated measures of association}
#' \item{method}{The method used to find these estimates}
#' @export
#' @seealso [fitLinModels()]
#' @inheritParams sbivarSingle
#' @inheritParams buildWeightMat
#'
#' @examples
#' n=1e2;m=8e1;p=3;k=4
#' ims = 6
#' Xl = lapply(selfName(seq_len(ims)), function(i){n = rpois(1, n)
#'  matrix(rnorm(n*p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
#' })
#' Yl = lapply(selfName(seq_len(ims)), function(i){m = rpois(1, m)
#'  matrix(rnorm(m*k), m, k, dimnames = list(NULL, paste0("Y", seq_len(k))))
#' })
#' Cxl = lapply(Xl, function(x){n = nrow(x)
#'     matrix(runif(n*2), n, 2, dimnames = list(NULL, c("x", "y")))
#' })
#' Eyl = lapply(Yl, function(y){m = nrow(y)
#'     matrix(runif(m*2), m, 2, dimnames = list(NULL, c("x", "y")))
#' })
#' estGAMs = sbivarMulti(Xl, Yl, Cxl, Eyl, method = "GAMs")
#' estMoran = sbivarMulti(Xl, Yl, Cxl, Eyl, method = "Moran")
#' estCorrelations = sbivarMulti(Xl, Yl, Cxl, Eyl, method = "Correlation")
sbivarMulti = function(Xl, Yl, Cxl, Eyl, families = list("X" = gaussian(), "Y" = gaussian()),
                       method = c("GAMs", "Correlation", "Moran's I"), mapToFinest = FALSE,
                        wo = c("distance", "nn"), numNN = 8, n_points_grid = 6e2){
    method = match.arg(method)
    wo = match.arg(wo)
    stopifnot(is.numeric(numNN), is.numeric(n_points_grid), is.character(wo),
              is.character(method), all(vapply(families, FUN.VALUE = TRUE, is, "family")))
    foo = checkInputMulti(Xl, Yl, Cxl, Eyl)
    jointCoordinates <- missing(Eyl)
    out = if(method == "GAMs"){
        wrapGAMsMulti(Xl, Yl, Cxl, Eyl, families = families,
                 n_points_grid = n_points_grid)
    }  else if(method == "Correlation"){
        wrapCorrelationsMulti(Xl, Yl, Cxl, Eyl, mapToFinest = mapToFinest,
                              jointCoordinates = jointCoordinates)
    } else if (method == "Moran's I"){
        wrapMoransIMulti(Xl, Yl, Cxl, Eyl, wo = wo, numNN = numNN)
    }
    return(list("estimates" = out, "method" = method))
}

