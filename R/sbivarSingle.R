#' Test for bivariate spatial association in a single image
#'
#' The tests can be applied to either disjoint or joint coordinate sets.
#'
#' @details If only Cx is supplied and X and Y have the same number of rows, a joint analysis is performed
#' If Cx and Ey are provided, and X and Y have the same number of rows, equality of Cx and Ey is checked.
#' If true, a joint analysis is run, with a warning.
#'
#' @param X,Y Matrices of omics measurements
#' @param Cx,Ey Corresponding coordinate matrices of dimension two
#' @param method A character string, indicating which method to apply
#' @param gpParams Parameters of the Gaussian processes, see details
#' @param GPmethod,Quants,numLscAlts,optControl,corStruct Passed onto \link{fitGP}
#' @param n_points_grid,families Passed onto \link{wrapGAMs}
#' @param wo,numNN,eta passed onto \link{buildWeightMat}
#' @param cutoof,width Cutoff and width of the variogram estimation, passed onto \link[gstat]{variogram}
#' @param model Variogram model, passed onto \link[gstat]{vgm}
#' @param verbose Should info on type of analysis be printed?
#'
#' @details Any normalization of the data should happen prior to calling this function.
#' For instance, count data or metabolome data are best scaled to relative values and log-normalized prior to fitting GPs.
#' For GAMs, usually no normalization is needed, as the non-gaussianity is taken care of by
#' the outcome distribution, offset and link functions. Currently, identity, inverse and log-link are implemented.
#'
#' @returns A list with at least the following components
#' \item{result}{A matrix which contains at least a p-values ("pVal") and a
#' Benjamini-Hochberg adjusted p-value ("pAdj"), sorted by increasing p-value.}
#' \item{families}{As provided}
#' \item{method}{As provided}
#' \item{multi}{FALSE, a flag for the type of analysis}
#' @importFrom stats p.adjust
#' @importFrom methods is
#' @importFrom nlme corGaus lmeControl
#' @importFrom BiocParallel bpparam
#' @note All methods use multithreading on the cluster provided using the BiocParallel package
#' @details gpParams must be a list of length 2 with names 'X' and 'Y', consisting of matrices
#' with rownames "mean", "nugget", "range" and "sigma", and column names as in X and Y.
#' This argument allows to pass parameters of the Gaussian processes estimated with other software
#' to perform the score test.
sbivarSingle = function(X, Y, Cx, Ey, method = c("Moran's I", "GAMs", "Modified t-test", "GPs"),
      n_points_grid = 6e2, families = list("X" = gaussian(), "Y" = gaussian()),
      GPmethod = c("REML", "ML"), wo = c("exp", "distance", "nn"), numNN = 8, model = "Gau",
      gpParams, Quants = c(0.005, 0.5), numLscAlts = 10, width = cutoff/15, eta = 0.05, cutoff = sqrt(2)/4,
      optControl = lmeControl(opt = "optim", maxIter = 5e2, msMaxIter = 5e2,
                              niterEM = 1e3, msMaxEval = 1e3),
      corStruct = corGaus(form = ~ x + y, nugget = TRUE, value = c(1, 0.25)), verbose = TRUE){
    stopifnot(is.numeric(n_points_grid), ncol(Cx) == 2, is.character(method), is.numeric(numNN),
              all(vapply(families, FUN.VALUE = TRUE, is, "family")), is.list(optControl),
             inherits(corStruct, "corStruct"), inherits(corStruct, "corGaus"),
             length(Quants)==2, is.numeric(Quants), is.logical(verbose))
    if(verbose){
        message("Starting sbivar analysis of a single image on ",
        bpparam()$workers, " computing cores")
    }
    n = nrow(X);m = nrow(Y);p = ncol(X);k=ncol(Y)
    X = giveValidNames(X);Y = giveValidNames(Y)
    method = match.arg(method)
    GPmethod = match.arg(GPmethod)
    wo = match.arg(wo)
    foo = checkInputSingle(X, Y, Cx, Ey)
    if(missing(Ey)){
        #Run a joint analysis
        message("Only one coordinate matrix Cx supplied, and dimensions of X and Y do match.
             Performing an analysis with joint coordinate sets.")
        Ey = Cx
    } else if(method=="Modified t-test"){
       stop("Two coordinate matrices supplied, whereas the modified t-test requires only one.
            If coordinates are shared, omit Ey. If coordinates are disjoint, consider using method = 'MoransI'.")
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
    if(!missing(gpParams)){
        if(!is.list(gpParams) || names(gpParams) != c("X", "Y") ||
           !all(vapply(gpParams, FUN.VALUE = TRUE, function(x){
               identical(sort(rownames(x), c("mean", "nugget", "range", "sigma")))
               })) || !all(colnames(X) %in% colnames(gpParams$X) ||
                           !all(colnames(Y) %in% colnames(gpParams$Y)))){
            stop("gpParams must be a list with names 'X' and 'Y' consisting of
            matrices with rownames 'mean', 'nugget', 'range', 'sigma',
                 and colnames the same as matrices X and Y!")
        }
    }
    colnames(Cx) = colnames(Ey) = c("x", "y")
    out = if(method=="Moran's I"){
        (moranRes <- wrapMoransI(X = X, Y = Y, Cx = Cx, Ey = Ey, wo = wo, numNN = numNN,
                    eta = eta, width = width, verbose = verbose, cutoff = cutoff, model = model))$out
    } else if(method == "GAMs"){
        wrapGAMs(X = X, Y = Y, Cx = Cx, Ey = Ey, families = families,
                 n_points_grid = n_points_grid, verbose = verbose)
    } else if(method == "GPs"){
        wrapGPs(X = X, Y = Y, Cx = Cx, Ey = Ey, gpParams = gpParams, Quants = Quants,
                GPmethod = GPmethod, corStruct = corStruct, optControl = optControl,
                numLscAlts = numLscAlts, verbose = verbose)
    } else if(method == "Modified t-test"){
        wrapModTtest(X = X, Y = Y, Cx = Cx, verbose = verbose)
    }
    out = cbind(out, "pAdj" = p.adjust(out[, "pVal"], method = "BH"))
    result = out[order(out[, "pVal"]),]
    list("result" = result, "families" = families, "method" = method,
         "multi" = FALSE, "maxIxy" =  if(method=="Moran's I"){moranRes$maxIxy})
}

