library(testthat)
library(sbivar)
library(BiocParallel)
n=8e1;m=1e2;p=4;k=3
X = matrix(rnorm(n*p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
Y = matrix(rnorm(m*k), m, k, dimnames = list(NULL, paste0("Y", seq_len(k))))
Cx = matrix(runif(n*2), n, 2)
Ey = matrix(runif(m*2), m, 2)
colnames(Cx) = colnames(Ey) = c("x", "y")
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
test_check("sbivar")
