#' @description extractResultsMulti() returns the results as matrix, including adjusted p-values, and sorted by p-value.
#'
#' @param method Multiplicity correction method passed onto p.adjust
#'
#' @importFrom stats p.adjust anova
#' @importFrom methods is
#' @return For extractResultsMulti() a list of matrices, all containing estimate, standard error,
#' p-value and adjusted p-value
#' @export
#' @rdname fitLinModels
#' @order 2
extractResultsMulti <- function(result, designDf, method = "BH") {
    categoricalVars = getDiscreteVars(designDf)
    if(result$method=="Moran's I")
        wParChar = switch(result$wo, "Gauss" = "eta_", "nn" = "numNN_")
    outout = loadBalanceBplapply(selfName(names(result$result)), function(nm){
        res = result$result[[nm]]
        if(result$method=="Moran's I")
            names(res) = paste0(wParChar, names(res))
        id <- vapply(res, FUN.VALUE = TRUE, function(x) {
            is(x, "lmerModLmerTest") || is(x, "lm")
        })
        out <- if (!any(id)) {
            #warning("No lmerModLmerTest or lm models found!")
            return(list("Intercept" = NULL, "fixedEffects" = NULL))
        } else {
            Summaries <- lapply(res[id], function(x) summary(x)$coef)
            ints <- t(vapply(Summaries, FUN.VALUE = double(3), function(x) {
                x["(Intercept)", c("Estimate", "Std. Error", "Pr(>|t|)")]
            }))
            colnames(ints) <- c("Estimate", "SE", "pVal")
            AnovaTabs <- lapply(res[id], anova)
            fixedVars = selfName(unique(unlist(lapply(res[id], function(x){
                all.vars(terms(x))[-1]
            }))))
            fixedOut <- lapply(fixedVars, function(Var) {
                if(discr <- Var %in% categoricalVars){
                    unVals <- unique(designDf[[Var]])
                    emptyCoef <- rep_len(NA, length(unVals))
                    names(emptyCoef) <- paste0(Var, unVals) # Prepare empty coefficient
                } else {
                    emptyCoef <- double(1)
                    names(emptyCoef) <- Var
                }
                pVal <- vapply(AnovaTabs, FUN.VALUE = double(1), function(x) x[Var, "Pr(>F)"])
                coefs <- lapply(Summaries, function(coefObj) {
                    # Prepare the empty coefficient vector with all levels present. If
                    # outcome is NA for all levels, the factor level gets dropped,
                    # causing problems downstream.
                    if(discr){
                        rn <- intersect(names(emptyCoef), rownames(coefObj))
                        emptyCoef[rn] <- coefObj[rn, "Estimate"]
                        emptyCoef
                    } else {
                        coefObj[Var, "Estimate"]
                    }
                })
                coefMat <- matrix(unlist(coefs), byrow = TRUE, nrow = length(pVal))
                colnames(coefMat) <- names(emptyCoef)
                if (ncol(coefMat) == 2) {
                    # If two levels, the two coefficients are each other's opposite
                    # in sum coding, and NA can be replaced by the negative. For more than two
                    # levels, it's more complicated as NA may indicate that the parameter was not estimated
                    oneNA <- which(rowSums(naMat <- is.na(coefMat)) == 1)
                    for (i in oneNA) {
                        coefMat[i, naMat[i, ]] <- -coefMat[i, !naMat[i, ]]
                    }
                }
                cbind(coefMat, pVal = pVal)
            })
            if(result$method=="Moran's I"){
                intRes = c(ints[, "Estimate"], "pVal" = CCT(ints[, "pVal"]))
                fixRes = lapply(fixedOut, function(x) {
                    coefs = c(x[, colnames(x)!="pVal"])
                    names(coefs) = makeNames(colnames(x)[colnames(x)!="pVal"], rownames(x))
                    c(coefs, "pVal" = CCT(x[, "pVal"]))
                })
            } else {
                intRes = ints[1,]
                fixRes = lapply(fixedOut, function(x) x[[1]])
            }
            c(list("Intercept" = intRes), fixRes)
        }
    })
    fixResOut = lapply(selfName(names(outout[[1]])), function(ii){
        tmpMat = t(vapply(outout, FUN.VALUE = outout[[1]][[ii]], function(x) x[[ii]]))
        tmpMat = tmpMat[order(tmpMat[, "pVal"]),]
        cbind(tmpMat, "pAdj" = p.adjust(tmpMat[, "pVal"], method = method))
    })
    return(c(list("result" = fixResOut), result[intersect(names(result),
        c("method", "families", "wo", "wParams", "multi", "assayX", "assayY"))]))
}
