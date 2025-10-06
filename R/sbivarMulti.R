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
#' \item{multiplicity}{"multi", a flag for the type of analysis}
#' @seealso \link{fitLinModels}
#' @note All methods use multithreading on the cluster provided using the BiocParallel package
#' @inheritParams sbivarSingle
sbivarMulti = function(Xl, Yl, Cxl, Eyl, families = list("X" = gaussian(), "Y" = gaussian()),
                       method = c("GAMs", "Correlation", "Moran's I"), mapToFinest = FALSE,
                        wo = c("distance", "nn"), numNN = 8, n_points_grid = 6e2, verbose = TRUE){
    method = match.arg(method)
    wo = match.arg(wo)
    stopifnot(is.numeric(numNN), is.numeric(n_points_grid), is.character(wo), is.logical(verbose),
              is.character(method), all(vapply(families, FUN.VALUE = TRUE, is, "family")))
    if(verbose){
        message("Performing sbivar analysis on ", length(Xl), " images")
    }
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
    return(list("estimates" = out, "method" = method, "multiplicity" = "multi",
                "families" = families))
}

