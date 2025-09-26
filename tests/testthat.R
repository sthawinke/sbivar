library(testthat)
library(sbivar)
library(BiocParallel)
n=8e1;m=1e2;p=4;k=3
X = matrix(rnorm(n*p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
Y = matrix(rnorm(m*k), m, k, dimnames = list(NULL, paste0("Y", seq_len(k))))
Cx = matrix(runif(n*2), n, 2)
Ey = matrix(runif(m*2), m, 2)
colnames(Cx) = colnames(Ey) = c("x", "y")
resModtTest = sbivarSingle(X, Y, Cx, Ey, method = "Modified")
#Multiple images
ims = 6
Xl = lapply(selfName(seq_len(ims)), function(i){n = rpois(1, n)
 matrix(rnorm(n*p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
})
Yl = lapply(selfName(seq_len(ims)), function(i){m = rpois(1, m)
 matrix(rnorm(m*k), m, k, dimnames = list(NULL, paste0("Y", seq_len(k))))
})
Cxl = lapply(Xl, function(x){n = nrow(x)
    matrix(runif(n*2), n, 2, dimnames = list(NULL, c("x", "y")))
})
Eyl = lapply(Yl, function(y){m = nrow(y)
    matrix(runif(m*2), m, 2, dimnames = list(NULL, c("x", "y")))
})
toyDesign = data.frame("covariate" = rnorm(ims),
                       "cofactor" = factor(rep(c(TRUE, FALSE), length.out = ims)),
                       "group" = rep(c("control", "treatment"), length.out = ims))
# Register the parallel backend
nCores <- 2
if(.Platform$OS.type == "unix"){
    #On unix-based systems, use MulticoreParam
    register(MulticoreParam(nCores))
} else {
    #On windows, use makeCluster
    library(doParallel)
    Clus = makeCluster(nCores)
    registerDoParallel(Clus)
    register(DoparParam(), default = TRUE)
}
#register(SerialParam()) # Switch on when mapping test coverage
estGAMs = sbivarMulti(Xl, Yl, Cxl, Eyl, method = "GAMs")
estMoran = sbivarMulti(Xl, Yl, Cxl, Eyl, method = "Moran")
estCorrelations = sbivarMulti(Xl, Yl, Cxl, Eyl, method = "Correlation")
multiFitGams = fitLinModels(estGAMs, design = toyDesign, Formula = out ~ covariate + cofactor + (1|group))
multiFitMoran = fitLinModels(estMoran, design = toyDesign, Formula = out ~ covariate + cofactor)
multiFit = fitLinModels(estCorrelations, design = toyDesign, Formula = out ~ covariate + (1|group))
#Extract the results
resGams = extractResultsMulti(multiFitGams, design = toyDesign)
test_check("sbivar")
