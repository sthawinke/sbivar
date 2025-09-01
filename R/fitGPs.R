fitGPs = function(x, xMat, method = c("REML", "ML", "gaussian", "poisson", "gamma", "negbin"),
                   training_iter = 100L, device = "cuda", offset, beta){#, startLengthScale = max(dist(xMat))/4
    method = match.arg(method)
    if(method == "gaussian"){
        pyFit = fit_spatial_gp(xMat, x, training_iter = training_iter, device = device)
        #Now built into python
        return(unlist(pyFit$Pars))
    } else if(method == "poisson"){
        pyFit = fit_poisson_gp(xMat, x, training_iter = training_iter,
                               device = device, log_offset = offset, beta = beta)
        return(unlist(pyFit$Pars))
    } else if(method == "negbin"){
        pyFit = fit_nb_gp(xMat, x, training_iter = training_iter,
                          device = device, log_offset = offset, beta = beta)
        return(unlist(pyFit$Pars))
    } else if(method == "gamma"){
        pyFit = fit_gamma_gp(xMat, x, training_iter = training_iter,
                             device = device, log_offset = offset, beta = beta)
        return(unlist(pyFit$Pars))
    } else {
        xModGls <- fitGLS(data.frame("out" = x, xMat), method = method)
        #xSig = corMatrix(xModGls$modelStruct$corStruct, corr = TRUE)
        xPars = c(coef(xModGls$modelStruct$corStruct, unconstrained = FALSE), sigma(xModGls))
        names(xPars) = c("range", "nugget", "sigma")
        c(xPars, "mean" = coef(xModGls)[1])#"Sigma" = xSig,
    }
}
fitManyGPs = function(mat, coord, nCores = 1, method = "gaussian",  ...){#startLengthScale = getStartLengthScale(coord),
    cns = colnames(mat);names(cns) = cns
    if(method %in% c("poisson", "gamma", "negbin")){
        offset = log(rowSums(mat))
        betas = log(colSums(mat)/sum(mat))
    }
    simplify2array(mclapply(mc.cores = nCores, cns, function(cn){
        prepGLS(mat[, cn], xMat = coord, method = method, offset = offset,
                beta = betas[cn], ...)
    }))
}
