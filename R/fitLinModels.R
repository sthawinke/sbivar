#' Fit linear models on measures calculated for replicated images, and extract the results
#' @description Given measures estimated from replicated images using \link{sbivarMulti},
#' fit linear models to determine significance.
#' To maintain interpretability of the intercept, continuous fixed effect variables are centered,
#' and sum coding is used for the categorical ones.
#'
#' @param result A result with the measures of bivariate spatial association,
#' from a call to the \link{sbivar} function with multiple images
#' @param designDf A design dataframe
#' @param Formula A formula for the linear model to be fitted, can contain random effects.
#' @param Control A control list for lmerTest::lmer
#' @param inverseWeigh A boolean, should estimates be inverse weighed by variance for GAMs and Moran's I?
#' @param scaleByMax A boolean, should Moran's I be scaled by maximum values before plugging into the linear model
#' @param verbose Should a message with number of linear models and cores be printed?
#' @returns For fitLinModels(), a list of linear models
#' @export
#' @details The left hand side of "Formula" can be provided or not,
#' but it will be overridden by "out" for downstream analysis
#'
#' @examples
#' #' #Multi-image analysis on Vicari data, using GAMs
#' data(Vicari)
#' #Subset to 5 images and 500 spots limit computing time
#' VicariMultiTest <- lapply(Vicari, function(x) lapply(x[1:5], function(y) y[1:500, ]))
#' VicariRes <- sbivar(VicariMultiTest$TranscriptOutcomes, VicariMultiTest$MetaboliteOutcomes,
#'     VicariMultiTest$TranscriptCoords, VicariMultiTest$MetaboliteCoords,
#'     normX = "rel", normY = "rel",
#'     method = "GAM"
#' )
#' mouse <- substr(names(Vicari$TranscriptOutcomes)[1:5], 1, 10)
#' designDf <- data.frame("mouse" = mouse) # The design matrix
#' multiGAMLmms <- fitLinModels(VicariRes, designDf, Formula = ~ (1 | mouse))
#' # Extract the results
#' resGAMsMulti <- extractResultsMulti(multiGAMLmms, designDf = designDf)
#' head(resGAMsMulti$result$Intercept)
#' @importFrom lmerTest lmer
#' @importFrom stats formula terms model.matrix
#' @importFrom lme4 lmerControl .makeCC isSingular lFormula mkReTrms
#' @importFrom reformulas findbars nobars
#' @importFrom methods is
#' @importFrom smoppix centerNumeric named.contr.sum loadBalanceBplapply
#' @importFrom BiocParallel bplapply bpparam
#' @seealso \link[lmerTest]{lmer}, \link[stats]{lm}, \link[sbivar]{sbivarMulti}, \link[stats]{p.adjust}
#' @order 1
fitLinModels <- function(
      result, designDf, Formula, verbose = TRUE, inverseWeigh = result$method != "Correlation", scaleByMax = TRUE,
      Control = lmerControl(
          check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL),
          check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
          check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)
      )
) {
    stopifnot(
        all(c("estimates", "method", "multi") %in% names(result)),
        is.logical(inverseWeigh), is.logical(scaleByMax), length(result$estimates) == nrow(designDf), is.data.frame(designDf),
        is.character(Formula) || is(Formula, "formula")
    )
    measures <- result$estimates
    if (!result$multi) {
        stop("Fitting linear models only makes sense for multi-image analyses!")
    }
    if (inverseWeigh && (result$method == "Moran's I") && !result$returnSEsMoransI) {
        stop("Inverse weighing only possible if the variances of Moran's I are included!
             Rerun sbivar() with estimateSEsMoransI=TRUE.")
    }
    if (inverseWeigh && (result$method == "Correlation")) {
        warning("Inverse weighing not available for method = 'Correlation'.
                Performing an unweighted analysis", immediate. = TRUE)
        inverseWeigh <- FALSE
    }
    namesFun <- switch(result$method,
        "Correlation" = names,
        rownames
    )
    Features <- selfName(unique(unlist(lapply(measures, function(x) namesFun(x$res))))) # All feature pairs present
    iter <- selfName(if (moran <- result$method == "Moran's I") {
        result$wParams
    } else {
        1
    })
    # Prepare arrays of outcomes and weights
    outArr <- array(0,
        dim = c(nrow(designDf), length(Features), length(iter)),
        dimnames = list(names(measures), Features, names(iter))
    )
    if (inverseWeigh) {
        weightsArr <- outArr
    }
    for (i in names(measures)) {
        outArr[i, namesFun(measures[[i]]$res), ] <- if (result$method == "Correlation") measures[[i]]$res else measures[[i]]$res[, seq_along(iter)]
        if (scaleByMax && (result$method == "Moran's I")) {
            outArr[i, namesFun(measures[[i]]$res), ] <- t(t(outArr[i, namesFun(measures[[i]]$res), ]) / measures[[i]]$maxIxy)
        }
        if (inverseWeigh && (result$method %in% c("Moran's I", "GAMs"))) {
            weightsArr[i, namesFun(measures[[i]]$res), ] <- 1 / measures[[i]]$res[, seq_along(iter) + length(iter)]^2
        }
    }
    # Prepare design matrices for linear model fitting
    Formula <- replaceLhs(formula(Formula))
    # Replace outcome variable by "out"
    MM <- length(findbars(Formula)) > 0
    fixedVars <- all.vars(nobars(Formula)[[3]])
    baseDf <- data.frame("out" = 0, centerNumeric(designDf))
    if (is.null(fixedVars)) {
        contrasts <- NULL
    } else {
        discreteVars <- selfName(intersect(getDiscreteVars(designDf), fixedVars))
        contrasts <- lapply(discreteVars, function(x) named.contr.sum)
    }
    if (MM) {
        ff <- lFormula(Formula,
            data = baseDf, contrasts = contrasts,
            na.action = na.omit
        )
        # Prepare random effects fitting
    }
    modMat <- model.matrix(nobars(Formula), baseDf, contrasts.arg = contrasts)
    # Fixed effects model matrix
    Assign <- attr(modMat, "assign")
    if (verbose) {
        message(
            "Fitting ", length(Features), if (MM) " mixed" else " fixed",
            " effects models on ", bpparam()$workers, " cores"
        )
    }
    models <- loadBalanceBplapply(Features, function(feat) {
        lapply(selfName(names(iter)), function(it) {
            baseDf$out <- outArr[, feat, it]
            out <- if (sum(id <- !is.na(baseDf$out)) < 3) {
                NULL
            } else {
                if (MM) {
                    ff$fr <- ff$fr[id, , drop = FALSE]
                    ff$X <- ff$X[id, , drop = FALSE]
                    attr(ff$X, "assign") <- Assign
                    ff$reTrms$Zt <- ff$reTrms$Zt[, id, drop = FALSE]
                }
                fitLinModel(
                    ff = ff, y = baseDf[id, "out"], Terms = terms(Formula),
                    modMat = modMat[id, , drop = FALSE],
                    weights = if (inverseWeigh) weightsArr[id, feat, it] / sum(weightsArr[id, feat, it]),
                    Control = Control, MM = MM, Assign = Assign
                )
            }
        })
    })
    return(c(
        list("result" = models),
        result[intersect(names(result), c("method", "families", "wo", "multi", "assayX", "assayY", "wParams"))]
    ))
}
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
                out <- lmerTest:::as_lmerModLT(out, devfun = devfun)
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
