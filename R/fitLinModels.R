#' Fit linear models on replicated image results
#' @description Given measures estimated from replicated images using \link{sbivarMulti},
#' fit linear models to determine significance.
#' To maintain interpretability of the intercept, continuous fixed effect variables are centered,
#' and sum coding is used for the categorical ones.
#'
#' @param measures A list of measures of bivariate spatial association
#' @param design A design dataframe
#' @param Formula A formula for the linear model to be fitted, can contain random effects
#' @param Control A control list for lmerTest::lmer
#' @returns A data frame containing, sorted by ascending p-value
#' @export
#'
#' @examples
#' example(sbivarMulti, "sbivar")
#' toyDesign = data.frame("covariate" = rnorm(ims),
#' "cofactor" = factor(rep(c(TRUE, FALSE), length.out = ims)),
#'  "group" = rep(c("control", "treatment"), length.out = ims))
#' multiFitGams = fitLinModels(estGAMs, design = toyDesign,
#' Formula = out ~ covariate + cofactor + (1|group))
#' multiFitMoran = fitLinModels(estMoran, design = toyDesign,
#' Formula = out ~ covariate + cofactor)
#' multiFit = fitLinModels(estCorrelations, design = toyDesign, Formula =
#' out ~ covariate + (1|group))
#' #Extract the results
#' resGams = extractResultsMulti(multiFitGams, design = toyDesign)
#' head(resGams$Intercept)
#' @importFrom lmerTest lmer
#' @importFrom stats formula terms model.matrix
#' @importFrom lme4 lmerControl .makeCC isSingular lFormula mkReTrms findbars nobars
#' @importFrom methods is
#' @importFrom smoppix centerNumeric named.contr.sum loadBalanceBplapply
#' @importFrom BiocParallel bplapply
#' @seealso \link[lmerTest]{lmer}, \link[stats]{lm}, \link{sbivarMulti}
fitLinModels = function(measures, design, Formula, Control = lmerControl(
    check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL),
    check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),
    check.conv.hess = .makeCC(action = "ignore", tol = 1e-06))){
    stopifnot(names(measures) == c("estimates", "method"))
    method = measures$method;measures = measures$estimates
    stopifnot(length(measures)==nrow(design), is.data.frame(design),
              is.character(Formula) || is(Formula, "formula"))
    withWeights <- (method %in% c("GAMs"))#, "Correlation"
    namesFun = if(withWeights) rownames else names
    Features = selfName(unique(unlist(lapply(measures, namesFun)))) # All feature pairs present
    #Prepare matrices of outcomes and weights
    outMat = matrix(0, nrow = nrow(design), ncol = length(Features),
                    dimnames = list(names(measures), Features))
    if(withWeights){
        weightsMat = outMat
        for(i in names(measures)){
            weightsMat[i,namesFun(measures[[i]])] = 1/measures[[i]][, "se"]^2
        }
    }
    for(i in names(measures)){
        outMat[i,namesFun(measures[[i]])] = if(withWeights) measures[[i]][, "est"] else measures[[i]]
    }
    #Prepare design matrices for linear model fitting
    Formula = formula(Formula)
    MM <- length(findbars(Formula)) > 0
    fixedVars = all.vars(nobars(Formula)[[3]])
    baseDf = data.frame("out" = 0, centerNumeric(design))
    if (is.null(fixedVars)) {
        contrasts <- NULL
    } else {
        discreteVars <- selfName(intersect(getDiscreteVars(design), fixedVars))
        contrasts <- lapply(discreteVars, function(x) named.contr.sum)
    }
    if(MM){
        ff <- lFormula(Formula, data = baseDf, contrasts = contrasts, na.action = na.omit)
        #Prepare random effects fitting
    }
    modMat <- model.matrix(nobars(Formula), baseDf, contrasts.arg = contrasts) #Fixed effects model matrix
    Assign <- attr(modMat, "assign")
    models <- loadBalanceBplapply(Features, function(feat) {
        baseDf$out = outMat[, feat]
        out <- if (sum(id <- !is.na(baseDf$out)) < 3) {
            NULL
        } else {
            if(MM){
                ff$fr <- ff$fr[id,, drop = FALSE];ff$X <- ff$X[id,, drop = FALSE]
                attr(ff$X, "assign") <- Assign
                ff$reTrms$Zt <- ff$reTrms$Zt[, id, drop = FALSE]
            }
            fitLinModel(ff = ff, y = baseDf[id, "out"], Terms = terms(Formula),
                        modMat = modMat[id,,drop = FALSE],
                        weights = if(withWeights) weightsMat[id, feat]/sum(weightsMat[id, feat]),
                        Control = Control, MM = MM, Assign = Assign)
        }
        return(out)
    })
    return(models)
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
#' @details Code is based on smoppix:::fitSingleLmmModel, but may diverge
#'
#' @returns A fitted lmer or lm model
#' @importFrom lme4 mkLmerDevfun optimizeLmer mkMerMod
#' @importFrom stats lm.wfit lm.fit
#' @importFrom smoppix lm_from_wfit
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
        Fit = if(is.null(weights)){
                lm.fit(y = y, x = modMat)
            } else {
                lm.wfit(y = y, x = modMat, w = weights)
            }
        mod <- lm_from_wfit(Fit, y = y, Assign = Assign, Terms = Terms)
    }
    return(mod)
}
