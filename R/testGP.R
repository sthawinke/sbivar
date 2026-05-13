#' Perform a score test on the bivariate spatial association in a Gaussian process.
#'
#' This function tests for the variance of a random effect, causing the covariance, to be zero.
#' It is a score test as developed by \insertCite{Zhang2003}{sbivar}, with the test statistic
#' having a scaled chi-square distribution.
#'
#' @inheritParams sbivarSingle
#' @param x,y outcome vectors
#' @param crossBlocks An \eqn{n \times m \times L} array of cross-blocks \eqn{C_l} from
#'   \code{\link{buildAltSigmas}}
#' @param distMat The complete distance matrix of Cx and Ey stacked. Only needed when
#'   \code{sx}, \code{sy}, \code{derivX}, or \code{derivY} are not supplied.
#' @param solXonly,solYonly Parameters of the Gaussian processes of x and y
#' @param sx,sy Inverses of the covariance matrices of x and y respectively.
#'   Computed from \code{distMat} and \code{solXonly}/\code{solYonly} when missing.
#' @param derivX,derivY \eqn{n \times n \times 3} / \eqn{m \times m \times 3} arrays of
#'   covariance-parameter derivative matrices. Computed when missing.
#' @inheritParams buildAltSigmas
#' @details Two tests are performed, one for positive and one for negative association,
#' and two times the smallest p-value to achieve a two-sided test. The sign indicates
#' which direction was most significant.
#'
#' @returns A vector of length 2: a p-value and an indicator of the sign:
#' +1 for positive association, -1 for negative
#' @importFrom stats pchisq
#' @importFrom abind abind
#' @importFrom Matrix bdiag
#' @references
#' \insertAllCited{}
testGP <- function(x, y, crossBlocks, solXonly, solYonly,
                   sx, sy, derivX, derivY, distMat) {
    n <- dim(crossBlocks)[1]
    m <- dim(crossBlocks)[2]
    idN <- seq_len(n)
    idM <- n + seq_len(m)

    # ---- Compute missing precomputed quantities from distMat ----
    if (missing(sx)) {
        sx <- base::solve(buildSigmaGp(solXonly, distMat = distMat[idN, idN]))
    }
    if (missing(sy)) {
        sy <- base::solve(buildSigmaGp(solYonly, distMat = distMat[idM, idM]))
    }
    if (missing(derivX)) {
        derivX <- buildDerivArray(solXonly, distMat[idN, idN], "X")
    }
    if (missing(derivY)) {
        derivY <- buildDerivArray(solYonly, distMat[idM, idM], "Y")
    }

    # ---- Data vectors ----
    diffVec <- c(x, y) - rep(c(solXonly["mean"], solYonly["mean"]), times = c(n, m))

    # ---- Block-wise projection matrix P (avoids forming the full invW) ----
    # invW = bdiag(sx, sy), regMat = cbind(1, c(rep(0,n), rep(1,m)))
    # P = invW - invW*regMat * solve(regMat^T*invW*regMat) * regMat^T*invW
    # Using invW block structure: invW*regMat[:,1] = [sx*1; sy*1], invW*regMat[:,2] = [0; sy*1]
    sx_1  <- as.vector(sx %*% rep(1, n))   # = rowSums(sx), n-vector
    sy_1  <- as.vector(sy %*% rep(1, m))   # = rowSums(sy), m-vector
    a     <- sum(sx_1)                      # 1^T sx 1
    b     <- sum(sy_1)                      # 1^T sy 1
    Binv  <- solve(matrix(c(a + b, b, b, b), 2L, 2L))
    A_top <- cbind(sx_1, 0)                 # n x 2  (invW*regMat top rows)
    A_bot <- cbind(sy_1, sy_1)             # m x 2  (invW*regMat bottom rows)
    AB_top <- A_top %*% Binv               # n x 2
    AB_bot <- A_bot %*% Binv               # m x 2
    # P blocks
    P_nn <- sx - tcrossprod(AB_top, A_top)  # n x n
    P_mm <- sy - tcrossprod(AB_bot, A_bot)  # m x m
    P_nm <- -tcrossprod(AB_top, A_bot)      # n x m
    P    <- rbind(cbind(P_nn, P_nm), cbind(t(P_nm), P_mm))  # (n+m) x (n+m)

    # ---- e = 0.5 * tr(P) ----
    e <- 0.5 * (sum(diag(P_nn)) + sum(diag(P_mm)))

    # ---- invW*diffVec using block structure (no full invW needed) ----
    invWDiffVec <- c(as.vector(sx %*% diffVec[idN]),
                     as.vector(sy %*% diffVec[idM]))

    # ---- Efficient Fisher information for nuisance parameters ----
    derivListP <- list(
        "X" = arrayMatProd(M = P_nn, A = derivX),
        "Y" = arrayMatProd(M = P_mm, A = derivY)
    )
    Ithetatheta <- 0.5 * bdiagn(
        arrayProdTr(derivListP$X, derivListP$X),
        arrayProdTr(derivListP$Y, derivListP$Y)
    )
    sitt <- try(solve(Ithetatheta), silent = TRUE)
    idItt <- colnames(sitt)
    if (inherits(sitt, "try-error")) {
        idItt <- grep("range", colnames(Ithetatheta), invert = TRUE)
        sitt  <- try(solve(Ithetatheta[idItt, idItt]), silent = TRUE)
    }
    if (inherits(sitt, "try-error")) {
        idItt <- grep("sigma", colnames(Ithetatheta))
        sitt  <- try(solve(Ithetatheta[idItt, idItt]), silent = TRUE)
    }

    # ---- Score-test quantities via C++ (block-structure exploitation) ----
    internals  <- scoreTestInternals_cpp(P, crossBlocks, derivX, derivY, invWDiffVec)
    Itautau    <- internals$Itautau
    Itautheta  <- internals$Itautheta
    colnames(Itautheta) <- c(dimnames(derivX)[[3]], dimnames(derivY)[[3]])
    UPos       <- internals$UPos
    UNeg       <- internals$UNeg

    ItautauTilde <- Itautau - rowSums(
        (Itautheta[, idItt, drop = FALSE] %*% sitt) * Itautheta[, idItt, drop = FALSE]
    )
    kappaEst <- ItautauTilde / (2 * e)
    nu       <- 2 * e^2 / ItautauTilde

    pVals <- cbind(
        pchisq(UPos / kappaEst, df = nu, lower.tail = FALSE),
        pchisq(UNeg / kappaEst, df = nu, lower.tail = FALSE)
    )
    ps    <- apply(pVals, 2, CCT)
    idMin <- which.min(ps)
    c("pVal" = min(2 * ps[idMin], 1), "sign" = if (idMin == 1L) 1L else -1L)
}

#' Build a 3-slice array of covariance-parameter derivatives
#'
#' @param fittedGP The fitted Gaussian process (vector of 4 parameters)
#' @param distMat Distance matrix for the relevant modality
#' @param suffix Character suffix to append to parameter names ("X" or "Y")
#' @return An \eqn{n \times n \times 3} array with 3rd-dim names
#'   \code{c("sigma<suffix>", "range<suffix>", "nugget<suffix>")}
#' @importFrom abind abind
buildDerivArray <- function(fittedGP, distMat, suffix) {
    arr <- abind(
        arrayDeriv(fittedGP, distMat, what = "sigma"),
        arrayDeriv(fittedGP, distMat, what = "range"),
        arrayDeriv(fittedGP, distMat, what = "nugget"),
        along = 3
    )
    dimnames(arr)[[3]] <- paste0(c("sigma", "range", "nugget"), suffix)
    arr
}

#' Construct the nxnxg/2 array of derivatives for a nxn matrix to the g/2 covariance matrix parameters.
#'
#' @param fittedGP The fitted Gaussian process (vector of 4 parameters)
#' @param what For which parameter is the derivative required?
#' @return The matrix of derivatives
arrayDeriv <- function(fittedGP, distMat, what) {
    switch(what,
        "sigma" = 2 * fittedGP["sigma"] * (exp(-(distMat / fittedGP["range"])^2) * (1 - fittedGP["nugget"]) +
            diag(fittedGP["nugget"], nrow(distMat))),
        "nugget" = {
            tmp <- -exp(-(distMat / fittedGP["range"])^2) * fittedGP["sigma"]^2
            diag(tmp) <- 0
            tmp
        },
        "range" = {
            tmp <- 2 * (1 - fittedGP["nugget"]) * fittedGP["sigma"]^2 *
                exp(-(distMat / fittedGP["range"])^2) * distMat^2 / fittedGP["range"]^3
            diag(tmp) <- 0
            tmp
        }
    )
}
