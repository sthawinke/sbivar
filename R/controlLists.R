#' Define default argument lists for the MoransISingle(), GPsingle() and GAMsSingle() functions
#'
#' @inheritParams sbivarSingle
moranControl = function(etas = c(5e-6, 2e-4, 2e-2), findMaxW = FALSE, returnSEsMoransI = TRUE,
                        variogramModels = c("Exp", "Lin"), width = cutoff / 15, cutoff = sqrt(2) / 3,
                        wo = c("Gauss", "nn"), numNNs = c(4, 8, 24)){
    list(etas = etas, findMaxW = findMaxW, returnSEsMoransI = returnSEsMoransI,
         variogramModels = variogramModels, width = width, cutoff = cutoff,
         wo = wo, numNNs = numNNs)
}
GPcontrol = function(GPmethod = c("REML", "ML"), gpParams, Quants = c(0.005, 0.5), numLscAlts = 5,
                     optControl = lmeControl(
                         opt = "optim", maxIter = 5e2, msMaxIter = 5e2,
                         niterEM = 1e3, msMaxEval = 1e3
                     ),
                     correlation = corGaus(form = ~ x + y, nugget = TRUE, value = c(0.9 * max(apply(Cx, 2, function(x) diff(range(x)))), 0.25))){
    list(GPmethod = GPmethod, gpParams = gpParams, Quants = Quants, numLscAlts = numLscAlts,
         optControl = optControl, correlation = correlation)

}
GAMcontrol = function(Gamm = FALSE, correlation = corGaus(form = ~ x + y, nugget = TRUE, value = c(0.9 * max(apply(Cx, 2, function(x) diff(range(x)))), 0.25))){
    list(Gamm = Gamm, correlation = correlation)
}
