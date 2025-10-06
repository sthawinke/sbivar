#' @description extractResultsMulti() returns the results as matrix, including adjusted p-values, and sorted by p-value.
#'
#' @param models The fitted linear models
#' @param method Multiplicity correction method passed onto p.adjust
#'
#' @importFrom stats p.adjust anova
#' @importFrom methods is
#' @return For extractResultsMulti() a list of matrices, all containing estimate, standard error,
#' p-value and adjusted p-value
#' @export
#' @rdname fitLinModels
#' @order 2
extractResultsMulti <- function(result, design, method = "BH") {
    id <- vapply(result$result, FUN.VALUE = TRUE, function(x) {
        is(x, "lmerModLmerTest") || is(x, "lm")
    })
    out <- if (!any(id)) {
        warning("No lmerModLmerTest or lm models found!")
        return(list("Intercept" = NULL, "fixedEffects" = NULL))
    } else {
        Summaries <- loadBalanceBplapply(result$result[id], function(x) summary(x)$coef)
        ints <- t(vapply(Summaries, FUN.VALUE = double(3), function(x) {
            x["(Intercept)", c("Estimate", "Std. Error", "Pr(>|t|)")]
        }))
        colnames(ints) <- c("Estimate", "SE", "pVal")
        intMat <- cbind(ints, pAdj = p.adjust(ints[, "pVal"], method = method))[order(ints[
            ,"pVal"]), ] # Order by p-value
        AnovaTabs <- loadBalanceBplapply(result$result[id], anova)
        fixedVars = selfName(unique(unlist(lapply(result$result[id], function(x){
            all.vars(terms(x))[-1]
        }))))
        categoricalVars = getDiscreteVars(design)
        fixedOut <- lapply(fixedVars, function(Var) {
            if(discr <- Var %in% categoricalVars){
                unVals <- unique(design[[Var]])
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
            cbind(coefMat, pVal = pVal, pAdj = p.adjust(pVal, method = method))[order(pVal), ]
        })
        c(list("Intercept" = intMat), fixedOut)
    }
    return(c(list("result" = out, "assayX" = result$assayX, "assayY" = result$assayY),
             result[c("method", "families", "multiplicity")]))
}
