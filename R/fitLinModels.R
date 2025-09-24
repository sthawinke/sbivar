#' Fit linear models on replicated image results
#' @description Given measures estimated from replicated images, fit linear models to determine significance.
#' The measures are usually estimated using \link{sbivarMulti}
#'
#' @param measures A list of measures of bivariate spatial association
#' @param design A design dataframe
#' @param Formula A formula for the linear model to be fitted, can contain random effects
#'
#' @returns A data frame containing, sorted by ascending p-value
#' @export
#'
#' @examples
#' example(sbivarMulti, "sbivar")
#' toyDesign = data.frame("covariate" = rnorm(im))
#' multiFitGams = fitLinModels(estGAMs, design = toyDesign,)
#' @importFrom lmerTest lmer
#' @importFrom stats formula terms model.matrix
#' @importFrom lme4 lmerControl .makeCC isSingular lFormula mkReTrms findbars
#' @importFrom methods is
#' @seealso \link[lmerTest]{lmer}, \link[stats]{lm}
fitLinModels = function(measures, design, Formula){
    stopifnot(length(measures)==nrow(design), is.matrix(design) || is.data.frame(design),
              is.character(Formula) || is(Formula, "formula"))
    Formula = formula(Formula)
    MM <- any(grepl("\\|", formChar <- characterFormula(Formula)))
    unPairs = unique(unlist(lapply(measures, names)))
    Control <- lmerControl(
        check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL),
        check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),
        check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)
    )
    #Matrix of outcomes
    outMat = matrix(0, nrow = nrow(design), ncol = length(unPairs),
                    dimnames = list(rownames(design), unPairs))
    for(i in names(measures)){
        measures[[i]][unPairs]
    }
    #Prepare design matrices
    design = centerNumeric(design)
    if(MM){
        ff <- lFormula(Formula, data = data.frame("out" = 0, design),
                       contrasts = contrasts, na.action = na.omit)
    }
    if (is.null(fixedVars)) {
        contrasts <- NULL
    } else {
        discreteVars <- selfName(intersect(getDiscreteVars(obj), fixedVars))
        contrasts <- lapply(discreteVars, function(x) named.contr.sum)
    }
    modMat <- model.matrix(nobars(Formula), baseDf, contrasts.arg = contrasts) #Fixed effects model matrix
    Assign <- attr(modMat, "assign")
    models <- loadBalanceBplapply(Features, function(gene) {
        df <- buildDataFrame(obj, gene = gene, pi = pi, pppDf = pppDf)
        out <- if (is.null(df) || sum(!is.na(df$pi)) < 3) {
            NULL
        } else {
            contrasts <- contrasts[!names(contrasts) %in% vapply(df, FUN.VALUE = TRUE, is.numeric)]
            fitLinModel(Formula, df, contrasts,
                       Control, MM = MM, Weight = df$weight)
        }
        return(out)
    })
}
#' Fit a linear model for an individual feature pair
#'
#' @param dff The dataframe
#' @param Control Control parameters
#' @param MM A boolean, should a mixed model be tried
#' @param Weight A weight variable
#' @inheritParams fitLinModels
#' @return A fitted model
#'
#' @details Code is based on smoppix:::fitPiModel, but may diverge
#'
#' @importFrom lmerTest lmer
#' @importFrom stats na.omit lm
#' @importFrom lme4 nobars
#' @seealso \link{fitLinModels}
fitLinModel <- function(Formula, dff, Control, MM, Weight = NULL) {
    environment(ff) <- environment() # Make sure formula sees the weights, see
    # https://stackoverflow.com/questions/61164404/call-to-weight-in-lm-within-function-doesnt-evaluate-properly
    if (MM) {
        mod <- try(lmerTest::lmer(ff, data = dff, na.action = na.omit,
                                  weights = Weight, control = Control), silent = TRUE)
    }
    # If no random variables or fit failed, go for linear model
    if (!MM || is(mod, "try-error")) {
        mod <- try(lm(nobars(ff), data = dff, na.action = na.omit,
                      weights = Weight), silent = TRUE)
    }
    return(mod)
}
