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
    contrasts <- contrasts[!names(contrasts) %in% vapply(df, FUN.VALUE = TRUE, is.numeric)]
    models <- loadBalanceBplapply(Features, function(gene) {
        df <- buildDataFrame(obj, gene = gene, pi = pi, pppDf = pppDf)
        out <- if (is.null(df) || sum(id <- !is.na(df$pi)) < 3) {
            NULL
        } else {
            if(MM){
                ff$fr <- ff$fr[id,, drop = FALSE];ff$X <- ff$X[id,, drop = FALSE]
                attr(ff$X, "assign") <- Assign
                ff$reTrms$Zt <- ff$reTrms$Zt[, id, drop = FALSE]
            }

            fitLinModel(ff = ff, y = mat[id, "pi"], Terms = terms(Formula), modMat = modMat[id,,drop = FALSE],
                        weights = mat[id, "weights"], Control = Control, MM = MM, Assign = Assign)
        }
        return(out)
    })
}
#' Fit a linear model for an individual feature pair
#'
#' @param ff The prepared frame
#' @param y outcome vector
#' @param weights weights vector
#' @param Assign,Terms Added to fitted fixed effects model
#' @param modMat Design matrix of the fixed effects model
#' @inheritParams fitLinModels
#' @return A fitted model
#'
#' @details Code is based on smoppix:::fitSingleLmmModel, but may diverge
#'
#' @returns A fitted lmer model
#' @importFrom lme4 mkLmerDevfun optimizeLmer mkMerMod
#' @importFrom stats lm.wfit
fitLinModel <-function(ff, y, Control, Terms, modMat, MM, Assign, weights = NULL) {
    if(MM){
        fr <- ff$fr                    # this is a data.frame (model frame)
        ## Use model-frame column names used by stats::model.frame
        fr$`(weights)` <- weights
        fr[["out"]] <- y               # replace response
        mod <- try({
            devfun <- mkLmerDevfun(fr, ff$X, ff$reTrms, control = Control)
            opt <- optimizeLmer(devfun, control = Control)
            out <- mkMerMod(rho = environment(devfun), opt = opt,
                            reTrms = ff$reTrms, fr = fr)
            out <- lmerTest:::as_lmerModLT(out, devfun = devfun)
        }, silent = TRUE)
    }
    # Switch to fixed effects model when fit failed
    if(!MM || inherits(mod, "try-error")){
        mod <- lm_from_wfit(lm.wfit(y = y, x = modMat, w = weights), y = y,
                            Assign = Assign, Terms = Terms)
    }
    return(mod)
}
