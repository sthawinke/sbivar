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
#' @param n_points_grid,families Passed onto \link{GAMsSingle}
#' @param wo,numNN passed onto \link{buildWeightMat}
#' @param etas A vector of decay parameters for the weight function, see details
#' @param cutoff,width Cutoff and width of the variogram estimation, passed onto \link[gstat]{vgm}
#' @param verbose Should info on type of analysis be printed?
#' @param findMaxW Is the maximum bivariate Moran's I needed?
#' @param pseudoCount A pseudocount added prior to log-normalization to avoid taking the log of zero
#' @param normX,normY Character vectors indicating normalization, "log" means log-normalization of relative abundances
#' @param variogramModels A character string, indicating the variogram model passed onto \link[gstat]{vgm}.
#' Currently, only "Exp" and "Lin" are implemented for computational reasons.
#' @param returnSEsMoransI A boolean, are standard errors of Moran's I to be returned?
#'
#' @details Any normalization of the data should happen prior to calling this function.
#' For instance, count data or metabolome data are best scaled to relative values and log-normalized prior to fitting GPs.
#' For GAMs, usually no normalization is needed, as the non-gaussianity is taken care of by
#' the outcome distribution, offset and link functions. Currently, identity, inverse and log-link are implemented.
#' If multiple decay parameters eta are supplied (etas is a vector), tests are performed for all weight matrices
#' and the resulting p-values combined using the Cauchy combination rule.
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
#' @seealso \link{MoransISingle}, \link{ModTtestSingle}, \link{GAMsSingle}, \link{GPsSingle}
sbivarSingle = function(X, Y, Cx, Ey, method = c("Moran's I", "GAMs", "Modified t-test", "GPs"),
      normX = c("none", "rel", "log"), normY = c("none", "rel", "log"),
      etas = c(2e-5, 2e-3, 2e-1), findMaxW = FALSE,
      families = list("X" = gaussian(), "Y" = gaussian()), n_points_grid = 6e2, verbose = TRUE,
      variogramModels = c("Exp", "Lin"), width = cutoff/15, cutoff = sqrt(2)/3,
      wo = c("Gauss", "nn"), numNN = c(4, 8, 24), pseudoCount = 1e-8, returnSEsMoransI = FALSE,
      GPmethod = c("REML", "ML"), gpParams, Quants = c(0.005, 0.5), numLscAlts = 5,
      optControl = lmeControl(opt = "optim", maxIter = 5e2, msMaxIter = 5e2,
                              niterEM = 1e3, msMaxEval = 1e3),
      corStruct = corGaus(form = ~ x + y, nugget = TRUE, value = c(1, 0.25))){
    stopifnot(is.numeric(n_points_grid), ncol(Cx) == 2, is.numeric(numNN), all(numNN>0),
              all(vapply(families, FUN.VALUE = TRUE, is, "family")), is.list(optControl),
             inherits(corStruct, "corStruct"), inherits(corStruct, "corGaus"), is.numeric(etas),
             length(Quants)==2, is.numeric(Quants), is.logical(verbose), is.logical(findMaxW))
    if(verbose){
        message("Starting sbivar analysis of a single image on ",
        bpparam()$workers, " computing cores")
    }
    X = addDimNames(X, "X");Y = addDimNames(Y, "Y")
    method = match.arg(method);variogramModels = match.arg(variogramModels, several.ok = TRUE)
    normX = match.arg(normX);normY = match.arg(normY)
    GPmethod = match.arg(GPmethod)
    wo = match.arg(wo)
    foo = checkInputSingle(X, Y, Cx, Ey)
    if(missing(Ey)){
        if(nrow(X)!=nrow(Y)){
            stop("Only one coordinate matrix supplied, but sample size of X and Y do not match!
                 Please supply the 'Ey' argument and consider 'Moran's I' as method.")
        }
        #Run a joint analysis
        if(verbose){
            message("Only one coordinate matrix Cx supplied, and dimensions of X and Y do match.
                 Performing an analysis with joint coordinate sets.")}
        Ey = Cx
    } else if(method=="Modified t-test"){
       stop("Two coordinate matrices supplied, whereas the modified t-test requires only one.
            If coordinates are shared, omit Ey. If coordinates are disjoint, consider using method = 'MoransI'.")
    }
    if(!identical(names(families), c("X", "Y"))){
        stop("Name families 'X' and 'Y' for unambiguous matching")
    }
    rownames(Cx) = rownames(X);rownames(Ey) = rownames(Y)
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
    if((normX == "log" || normY == "log") && method == "GAMs"){
        warning("Normalizing data is not recommended for GAMs!
                try accounting for non-normality through the 'families' argument.")
    }
    colnames(Cx) = colnames(Ey) = c("x", "y")
    X = normMat(X, normX, pseudoCount);Cx = Cx[rownames(X),]
    Y = normMat(Y, normY, pseudoCount);Ey = Ey[rownames(Y),]
    out = if(method=="Moran's I"){
        (moranRes <- MoransISingle(X = X, Y = Y, Cx = Cx, Ey = Ey, wo = wo, numNN = numNN,
            variogramModels = variogramModels, etas = selfName(etas), width = width,
            returnSEsMoransI = returnSEsMoransI, verbose = verbose, cutoff = cutoff, findMaxW = findMaxW))$res
    } else if(method == "GAMs"){
        GAMsSingle(X = X, Y = Y, Cx = Cx, Ey = Ey, families = families,
                 n_points_grid = n_points_grid, verbose = verbose)
    } else if(method == "GPs"){
        GPsSingle(X = X, Y = Y, Cx = Cx, Ey = Ey, gpParams = gpParams, Quants = Quants,
                GPmethod = GPmethod, corStruct = corStruct, optControl = optControl,
                numLscAlts = numLscAlts, verbose = verbose)
    } else if(method == "Modified t-test"){
        sharedNames = intersect(rownames(X), rownames(Y))
        ModTtestSingle(X = X[sharedNames,], Y = Y[sharedNames,],
                     Cx = Cx[sharedNames,], verbose = verbose)
    }
    out = cbind(out, "pAdj" = p.adjust(out[, "pVal"], method = "BH"))
    result = out[order(out[, "pVal"]),]
    lis = list("result" = result, "method" = method,
         "multi" = FALSE, "normX" = normX, "normY" = normY)
    if(method=="Moran's I")
        lis$maxIxy = moranRes$maxIxy
    if(method=="GAMs")
        lis$families = families
    return(lis)
}

