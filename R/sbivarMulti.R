#' Apply one or more tests for bivariate spatial association
#'
#' This function calculates measures of spatial association for every image.
#' The resulting estimates can then be analysed further using the \link{fitLinModels} function.
#'
#' @param Xl,Yl Lists of matrices of omics measurements
#' @param Cxl,Eyl Lists of corresponding coordinate matrices of dimension two
#' @param method A character string, indicating which method to apply
#' @param families,n_points_grid Passed onto \link{wrapGAMsMulti}
#' @param wo,numNN Passed onto \link{buildWeightMat}
#' @returns A list containing
#' \item{estimates}{The estimated measures of association}
#' \item{method}{The method used to find these estimates}
#' \item{multi}{TRUE, a flag for the type of analysis}
#' @seealso \link{fitLinModels}
#' @note All methods use multithreading on the cluster provided using the BiocParallel package
#' @inheritParams sbivarSingle
#' @inheritParams buildWeightMat
#' @importFrom BiocParallel bpparam
sbivarMulti = function(Xl, Yl, Cxl, Eyl, families = list("X" = gaussian(), "Y" = gaussian()),
                       method = c("Moran's I", "GAMs", "Correlation"), eta = 0.025, normX = c("none", "log"),
                       normY = c("none", "log"),
                        wo = c("exp", "nn"), numNN = 8, n_points_grid = 6e2, verbose = TRUE){
    method = match.arg(method)
    wo = match.arg(wo)
    stopifnot(is.numeric(numNN), is.numeric(n_points_grid), is.character(wo), is.logical(verbose),
              is.character(method), all(vapply(families, FUN.VALUE = TRUE, is, "family")))
    if(verbose){
        message("Starting sbivar analysis of ", length(Xl), " images on ",
                bpparam()$workers, " computing cores")
    }
    Xl = lapply(Xl, addDimNames, "X");Yl = lapply(Yl, addDimNames, "Y")
    Cxl = mapply(Cxl, Xl, SIMPLIFY = FALSE, FUN = tmpFun <- function(cx, x) {colnames(cx) = c("x", "y");rownames(cx) = rownames(x);cx})
    Eyl = mapply(Eyl, Yl, SIMPLIFY = FALSE, FUN = tmpFun)
    normX = match.arg(normX);normY = match.arg(normY)
    foo = checkInputMulti(Xl, Yl, Cxl, Eyl)
    if(normX=="log"){
        Xl = lapply(Xl, logNorm)
        Cxl = lapply(selfName(names(Cxl)), function(i){Cxl[[i]][rownames(Xl[[i]]),]})
    }
    if(normY=="log"){
        Yl = lapply(Yl, logNorm)
        Eyl = lapply(selfName(names(Eyl)), function(i){Eyl[[i]][rownames(Yl[[i]]),]})
    }
    out = if(method == "GAMs"){
        wrapGAMsMulti(Xl, Yl, Cxl, Eyl, families = families,
                 n_points_grid = n_points_grid, verbose = verbose)
    }  else if(method == "Correlation"){
        wrapCorrelationsMulti(Xl, Yl, verbose = verbose)
    } else if (method == "Moran's I"){
        wrapMoransIMulti(Xl, Yl, Cxl, Eyl, wo = wo, numNN = numNN, verbose = verbose, eta = eta)
    }
    return(list("estimates" = out, "method" = method, "multi" = TRUE,
                "families" = families, "normX" = normX, "normY" = normY))
}

