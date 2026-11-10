#' Fit a linear model for an individual feature pair
#'
#' @param ff The prepared frame
#' @param y outcome vector
#' @param weights weights vector
#' @param Assign,Terms Added to fitted fixed effects model
#' @param modMat Design matrix of the fixed effects model
#' @param MM a Boolean, should a mixed model fit be attempted?
#' @inheritParams fitLinModels
#' @return A fitted model
#'
#' @details The code is based on \link[smoppix]{fitSingleLmmModel}, but may diverge
#'
#' @returns A fitted lmer or lm model
#' @importFrom lme4 mkLmerDevfun optimizeLmer mkMerMod
#' @importFrom lmerTest as_lmerModLmerTest
#' @importFrom stats lm.wfit lm.fit
#' @importFrom smoppix lm_from_wfit
fitLinModel <- function(ff, y, Control, Terms, modMat, MM, Assign, weights = NULL) {
  if (MM) {
    fr <- ff$fr # this is a data.frame (model frame)
    ## Use model-frame column names used by stats::model.frame
    fr$`(weights)` <- weights
    fr[["out"]] <- y # replace response
    mod <- try(
      {
        devfun <- mkLmerDevfun(fr, ff$X, ff$reTrms, control = Control)
        opt <- optimizeLmer(devfun, control = Control)
        out <- mkMerMod(
          rho = environment(devfun), opt = opt,
          reTrms = ff$reTrms, fr = fr
        )
        out <- as_lmerModLmerTest(out)
      },
      silent = TRUE
    )
  }
  # Switch to fixed effects model when fit failed
  if (!MM || inherits(mod, "try-error")) {
    Fit <- if (is.null(weights)) {
      lm.fit(y = y, x = modMat)
    } else {
      lm.wfit(y = y, x = modMat, w = weights)
    }
    mod <- lm_from_wfit(Fit, y = y, Assign = Assign, Terms = Terms)
  }
  return(mod)
}
