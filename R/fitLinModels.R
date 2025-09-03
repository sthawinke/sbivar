#' Fit linear models on replicated image results
#' @description Given measures estimated from replicated images, fit linear models to determine significance.
#' The measures are usually estimated using \link{sbivarMulti}
#'
#' @param measures A list of measures of bivairate spatial association
#' @param design A design dataframe
#' @param formula A formula for the linear model to be fitted, can contain random effects
#'
#' @returns A data frame containing, sorted by ascending p-value
#' @export
#'
#' @examples
#' example(sbivarMulti, "sbivar")
#' toyDesign = data.frame("covariate" = rnorm(im))
#' multiFitGams = fitLinModels(estGAMs, design = toyDesign,)
#' @importFrom stats formula
#' @seealso \link[lmerTest]{lmer}, \link[stats]{lm}
fitLinModels = function(measures, design, Formula){
    stopifnot(length(measures)==nrow(design), is.matrix(design) || is.data.frame(design),
              is.character(Formula) || is(Formula, "formula"))
    Formula = formula(Formula)
    unPairs = unique(unlist(lapply(measures, names)))
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
#' @seealso \link{fitLMMsSingle}
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
