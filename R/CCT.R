#' An analytical p-value combination method using the Cauchy distribution.
#'
#' The \code{CCT} function takes in a numeric vector of p-values, a numeric
#' vector of non-negative weights, and return the aggregated p-value using Cauchy combination rule
#' The code was taken from the xihaoli/STAAR github repo, and adapted
#' @param pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @param weights a numeric vector of non-negative weights. If \code{NULL}, the
#' equal weights are assumed.
#' @return The aggregated p-value
#' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
#' with analytic p-value calculation under arbitrary dependency structures.
#' \emph{Journal of the American Statistical Association}, \emph{115}(529), 393-402.
#' (\href{https://doi.org/10.1080/01621459.2018.1554485}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(3), 410-421.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.01.002}{pub})
#' @importFrom stats pcauchy
CCT <- function(pvals, weights=NULL){
    #### check if all p-values are between 0 and 1
    if(any(pvals<0) || any(pvals>1)){
        stop("All p-values must be between 0 and 1!")
    }

    #### check if there are p-values that are either exactly 0 or 1.
    is.zero <- any(pvals==0)
    is.one <- any(pvals==1)
    if(is.zero && is.one){
        return(1)
        #stop("Cannot have both 0 and 1 p-values!")
    } else if(is.zero){
        return(0)
    } else if(is.one){
        warning("There are p-values that are exactly 1!")
        return(1)
    }
    #### check the validity of weights (default: equal weights) and standardize them.
    if(is.null(weights)){
        weights <- rep(1/length(pvals),length(pvals))
    }else if(length(weights)!=length(pvals)){
        stop("The length of weights should be the same as that of the p-values!")
    }else if(any(weights < 0)){
        stop("All the weights must be positive!")
    }else{
        weights <- weights/sum(weights)
    }
    #### check if there are very small non-zero p-values
    is.small <- (pvals < 1e-16)
    if (!any(is.small)){
        cct.stat <- sum(weights*tan((0.5-pvals)*pi))
    }else{
        cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
        cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
    }

    #### check if the test statistic is very large.
    if(cct.stat > 1e+15){
        pval <- (1/cct.stat)/pi
    }else{
        pval <- pcauchy(cct.stat, lower.tail = FALSE)
    }
    return(pval)
}
