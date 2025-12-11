#' An analytical p-value combination method using the Cauchy distribution.
#'
#' The \code{CCT} function takes in a numeric vector of p-values,
#' and returns the aggregated p-value using Cauchy combination rule
#' The code was taken from the xihaoli/STAAR github repo, and adapted.
#' @param pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @return The aggregated p-value
#' \insertAllCited{}
#' @importFrom stats pcauchy
CCT <- function(pvals){
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
    #### check if there are very small non-zero p-values
    is.small <- (pvals < 1e-16)
    if (!any(is.small)){
        cct.stat <- mean(tan((0.5-pvals)*pi))
    }else{
        cct.stat <- sum(1/(pi*pvals[is.small]))
        cct.stat <- (cct.stat + sum(tan((0.5-pvals[!is.small])*pi)))/length(pvals)
    }
    #### check if the test statistic is very large.
    pval <-  if(cct.stat > 1e+15){
        (1/cct.stat)/pi
    } else {
        pcauchy(cct.stat, lower.tail = FALSE)
    }
    return(pval)
}
