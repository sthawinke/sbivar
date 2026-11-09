#' @importFrom utils head
#' @rdname fitLinModels
setMethod(
    "fitLinModels",
    "sbivarResultsMulti",
    function(result, designDf, Formula, verbose = TRUE, inverseWeigh = FALSE, scaleByMax = TRUE,
    Control = lmerControl(
        check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL),
        check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
        check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)
    )) {
        if (missing(designDf) && length(ncol(attr(terms(Formula), "factors"))) == 0) {
            designDf <- data.frame("foo" = seq_along(result@result))
        } # Allow for intercept only models
        stopifnot(
            is.logical(inverseWeigh), is.logical(scaleByMax), length(result@result) == nrow(designDf),
            is.data.frame(designDf), is.character(Formula) || is(Formula, "formula")
        )
        if (!result@multi) {
            stop("Fitting linear models only makes sense for multi-image analyses!")
        }
        if (inverseWeigh && (result@method == "Moran's I") && !result@returnSEsMoransI) {
            stop(strwrap("Inverse weighing only possible if the variances of Moran's I are included!
             Rerun sbivar() with returnSEsMoransI=TRUE."))
        }
        if (inverseWeigh && (result@method == "Correlation")) {
            warning(strwrap("Inverse weighing not available for method = 'Correlation'.
                Performing an unweighted analysis"), immediate. = TRUE)
            inverseWeigh <- FALSE
        }
        namesFun <- switch(result@method,
            "Correlation" = names,
            rownames
        )
        Features <- selfName(unique(unlist(lapply(result@result, function(x) namesFun(x$res))))) # All feature pairs present
        iter <- selfName(if (moran <- result@method == "Moran's I") {
            result@wParams
        } else {
            1
        })
        # Prepare arrays of outcomes and weights
        outArr <- array(NA,
            dim = c(nrow(designDf), length(Features), length(iter)),
            dimnames = list(names(result@result), Features, names(iter))
        )
        if (inverseWeigh) {
            weightsArr <- outArr
        }
        for (i in names(result@result)) {
            outArr[i, namesFun(result@result[[i]]$res), ] <- if (result@method == "Correlation") result@result[[i]]$res else result@result[[i]]$res[, seq_along(iter)]
            if (scaleByMax && (result@method == "Moran's I")) {
                # Scale estimates
                outArr[i, namesFun(result@result[[i]]$res), ] <- t(t(outArr[i, namesFun(result@result[[i]]$res), ]) / result@result[[i]]$maxIxy)
            }
            if (inverseWeigh && (result@method %in% c("Moran's I", "GAMs"))) {
                weightsArr[i, namesFun(result@result[[i]]$res), ] <- 1 / result@result[[i]]$res[, seq_along(iter) + length(iter)]^2
                if (scaleByMax && (result@method == "Moran's I")) {
                    # Scale variances too
                    weightsArr[i, namesFun(result@result[[i]]$res), ] <- t(t(weightsArr[i, namesFun(result@result[[i]]$res), ]) * result@result[[i]]$maxIxy^2)
                }
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
                " effects models on ", bpworkers(bpparam()), " cores"
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
            result[intersect(names(result), c(
                "method", "families", "wo", "multi", "assayX", "assayY",
                "wParams", "normX", "normY", "estimateSEsMoransI", "Gamm", "correlation"
            ))]
        ))
    }
)
