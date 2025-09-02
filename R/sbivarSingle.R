#' Apply one or more tests for bivariate spatial association
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
#' @param families A vector of length 2 giving outcome values.
#' @param n_points_grid The number of points in the new grid for the GAMs to be
#' evaluated on.
#' @param mapToFinest A boolean, should the one-to-one mapping for modified t-test
#' occur to the dataset with the best resolution?
#' @inheritParams wrapGPs
#' @inheritParams fitGPs
#'
#' @details Any normalization of the data should happen prior to calling this function.
#' For instance, count data or metabolome data are best scaled to relative values and log-normalized prior to fitting GPs.
#'
#' @returns A matrix which contains at least a p-values ("pVal") and a Benjamini-Hochberg adjusted p-value ("pAdj"),
#' sorted by increasing p-value.
#' @export
#' @importFrom stats p.adjust
#' @importFrom methods is
#' @importFrom nlme corGaus lmeControl
#' @examples
#' n=1e2;m=2e2;p=10;k=5
#' X = matrix(rnorm(n*p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
#' Y = matrix(rnorm(m*k), m, k, dimnames = list(NULL, paste0("Y", seq_len(k))))
#' Cx = matrix(runif(n*2), n, 2)
#' Ey = matrix(runif(m*2), m, 2)
#' colnames(Cx) = colnames(Ey) = c("x", "y")
#' resGAMs = sbivarSingle(X, Y, Cx, Ey, method = "GAMs")
#' resModtTest = sbivarSingle(X, Y, Cx, Ey, method = "Modified")
#' resModtTestJoint = sbivarSingle(X, Y[seq_len(nrow(X)),], Cx, method = "Modified")
#' resModtGPs = sbivarSingle(X, Y, Cx, Ey, method = "GPs")
sbivarSingle = function(X, Y, Cx, Ey, method = c("GAMs", "Modified t-test", "GPs"),
                  n_points_grid = 5e2, mapToFinest = FALSE, families = list("X" = gaussian(), "Y" = gaussian()),
                  GPmethod = c("REML", "ML", "gpytorch"), device = c("cpu", "cuda"),
                  training_iter = 100L, gpParams, Quants = c(0.005, 0.5), numLscAlts = 10,
                  optControl = lmeControl(opt = "optim", maxIter = 5e2,
                                                msMaxIter = 5e2, niterEM = 1e3,
                                                msMaxEval = 1e3),
                  corStruct = corGaus(form = ~ x + y, nugget = TRUE, value = c(1, 0.25))){
    stopifnot(is.numeric(n_points_grid), ncol(Cx) == 2, is.character(method),
              all(vapply(families, FUN.VALUE = TRUE, is, "family")), is.list(optControl),
              is.numeric(training_iter), inherits(corStruct, "corStruct"),
              inherits(corStruct, "corGaus"), length(Quants)==2, is.numeric(Quants))
    n = nrow(X);m = nrow(Y);p = ncol(X);k=ncol(Y)
    method = match.arg(method)
    GPmethod = match.arg(GPmethod)
    device = match.arg(device)
    #Check if only Cx supplied, or Cx and Ey are identical
    if(missing(Ey)){
        if(n!=m){
            stop("Only one coordinate matrix Cx supplied, and dimensions of X and Y do not match.
                 Please provide the coordinates of Y too through the Ey argument.")
        } else {#Run a joint analysis
            message("Only one coordinate matrix Cx supplied, and dimensions of X and Y do match.
                 Performing a joint analysis.")
            Ey = Cx
            jointCoordinates = TRUE
        }
    } else if(identical(Cx, Ey)){
        message("Identical coordinate matrices supplied, performing a joint analysis")
        jointCoordinates = TRUE
    } else {
        jointCoordinates = FALSE
    }
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
    out = if(method == "GAMs"){
        wrapGAMs(X = X, Y = Y, Cx = Cx, Ey = Ey, families = families,
                 n_points_grid = n_points_grid)
    } else if(method == "GPs"){
        wrapGPs(X = X, Y = Y, Cx = Cx, Ey = Ey, gpParams = gpParams, Quants = Quants,
                GPmethod = GPmethod, device  = device, training_iter = training_iter,
                corStruct = corStruct, optControl = optControl, numLscAlts = numLscAlts)
    } else if(method == "Modified t-test"){
        wrapModTtest(X = X, Y = Y, Cx = Cx, Ey = Ey, mapToFinest = mapToFinest,
                     jointCoordinates = jointCoordinates)
    }
    out = cbind(out, "pAdj" = p.adjust(out[, "pVal"], method = "BH"))
    out[order(out[, "pVal"]),]
}

