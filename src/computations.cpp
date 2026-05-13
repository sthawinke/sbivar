#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Evaluate Exp or Lin variogram model on a vector of distances
//'
//' Faster C++ replacement for the R \code{evalVariogram} wrapper.
//' Implements the Exponential model (\code{psill * exp(-d / range)}) and
//' the Linear model (\code{psill * (1 - d / range)} for \code{d < range},
//' else 0).
//'
//' @param distVec Numeric vector of pairwise distances
//' @param psill   Partial sill of the structured spatial component
//' @param range_  Range parameter
//' @param modelExp Logical: \code{TRUE} for Exponential, \code{FALSE} for Linear
//' @return Numeric vector of covariance values, same length as \code{distVec}
//' @keywords internal
// [[Rcpp::export]]
arma::vec evalVariogramCpp(
    const arma::vec& distVec,
    double psill,
    double range_,
    bool modelExp
) {
    arma::uword n = distVec.n_elem;
    arma::vec out(n, arma::fill::zeros);
    if (modelExp) {
        out = psill * arma::exp(-distVec / range_);
    } else {
        // Linear (spherical-like): psill * (1 - d/range) for d < range, else 0
        for (arma::uword i = 0; i < n; i++) {
            if (distVec[i] < range_) {
                out[i] = psill * (1.0 - distVec[i] / range_);
            }
        }
    }
    return out;
}

//' Batch-compute lower triangles and traces of \eqn{W_i^T \Sigma_X W_i}
//'
//' For each weight-matrix slice \eqn{W_i} of the \eqn{n \times m \times numWs}
//' array \code{W}, constructs the spatial covariance matrix \eqn{\Sigma_X} from
//' pre-evaluated lower-triangle covariance values and computes
//' \eqn{W_i^T \Sigma_X W_i} (\eqn{m \times m}), returning its lower triangle and
//' trace.
//'
//' The lower-triangle iteration order (column-major) matches R's
//' \code{which(lower.tri(diag(n)))} so that \code{vgVals} can be passed directly
//' from \code{evalVariogramCpp(distX, ...)}, where \code{distX} is
//' \code{as.vector(stats::dist(Cx))}.
//'
//' @param vgVals Covariance values for the strictly lower triangle of
//'   \eqn{\Sigma_X}, in column-major order (length \eqn{n(n-1)/2}).
//'   The diagonal is set to 1 (unit variance after standardisation).
//' @param W \eqn{n \times m \times numWs} array of weight matrices
//' @return A list with
//'   \describe{
//'     \item{sigXws}{\eqn{mm2 \times numWs} matrix, where \eqn{mm2 = m(m-1)/2};
//'       column \eqn{i} holds the lower-triangle entries of \eqn{W_i^T \Sigma_X W_i}
//'       in column-major order.}
//'     \item{traces}{Length-numWs vector of \eqn{tr(W_i^T \Sigma_X W_i)}.}
//'   }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List computeSigXws(
    const arma::vec& vgVals,
    const arma::cube& W
) {
    int n   = W.n_rows;
    int m   = W.n_cols;
    int nWs = W.n_slices;
    int mm2 = m * (m - 1) / 2;

    // Build symmetric Sigma_X: diagonal = 1, off-diagonal from vgVals.
    // Column-major lower-triangle fill matches which(lower.tri(diag(n))) in R.
    arma::mat SigmaX(n, n, arma::fill::eye);
    {
        arma::uword idx = 0;
        for (int j = 0; j < n - 1; j++) {
            for (int i = j + 1; i < n; i++) {
                SigmaX(i, j) = vgVals[idx];
                SigmaX(j, i) = vgVals[idx];
                idx++;
            }
        }
    }

    arma::mat sigXws(mm2, nWs);
    // Use Rcpp::NumericVector so traces comes back as a plain R vector, not a matrix
    Rcpp::NumericVector traces(nWs);

    // Process one W slice at a time, avoiding the n x (m * nWs) intermediate SigW.
    // Peak additional memory is one n x m scratch matrix (SWi) instead of n x (m * nWs).
    for (int wi = 0; wi < nWs; wi++) {
        const arma::mat Wi(const_cast<double*>(W.slice_memptr(wi)), n, m, false, true);
        arma::mat SWi = SigmaX * Wi;   // n x m  — reused each iteration
        arma::mat tmp = Wi.t() * SWi;  // m x m

        // Extract strictly lower triangle in column-major order,
        // matching which(lower.tri(diag(m))) in R
        arma::uword k = 0;
        for (int j = 0; j < m - 1; j++) {
            for (int i = j + 1; i < m; i++) {
                sigXws(k++, wi) = tmp(i, j);
            }
        }
        traces[wi] = arma::trace(tmp);
    }

    return Rcpp::List::create(
        Rcpp::Named("sigXws") = sigXws,
        Rcpp::Named("traces") = traces
    );
}

//' Compute Itautau, Itautheta and score statistics for the GP score test
//'
//' For each length-scale slice \eqn{l}, exploits the block structure of
//' \eqn{\Sigma_{alt,l} = I + \begin{bmatrix}0 & C_l \\ C_l^T & 0\end{bmatrix}}
//' to compute \eqn{P \Sigma_{alt,l}} block-wise without materialising the full
//' \eqn{(n+m)^2 \times L} intermediate array.
//'
//' @param P \eqn{(n+m) \times (n+m)} projection matrix
//' @param crossBlocks \eqn{n \times m \times L} array of cross-blocks \eqn{C_l}
//'   (the off-diagonal blocks of each \eqn{\Sigma_{alt,l}})
//' @param derivX \eqn{n \times n \times 3} covariance-parameter derivative arrays for X
//' @param derivY \eqn{m \times m \times 3} covariance-parameter derivative arrays for Y
//' @param vecPos \eqn{(n+m)} score vector (\eqn{\Omega^{-1}(z - \mu)} for the positive direction)
//' @return A list with
//'   \describe{
//'     \item{Itautau}{Length-\eqn{L} vector of \eqn{0.5\,\|P\Sigma_{alt,l}\|_F^2}}
//'     \item{Itautheta}{\eqn{L \times 6} matrix of trace cross-products with derivative matrices}
//'     \item{UPos}{Length-\eqn{L} vector of score statistics for the positive direction}
//'     \item{UNeg}{Length-\eqn{L} vector of score statistics for the negative direction}
//'   }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List scoreTestInternals_cpp(
    const arma::mat&  P,
    const arma::cube& crossBlocks,
    const arma::cube& derivX,
    const arma::cube& derivY,
    const arma::vec&  vecPos
) {
    const arma::uword n  = crossBlocks.n_rows;
    const arma::uword m  = crossBlocks.n_cols;
    const arma::uword L  = crossBlocks.n_slices;
    const arma::uword nm = n + m;

    // Extract P blocks once — reused across all L slices
    const arma::mat P_nn = P.submat(0, 0,  n-1,  n-1);   // n x n
    const arma::mat P_nm = P.submat(0, n,  n-1, nm-1);   // n x m
    const arma::mat P_mm = P.submat(n, n, nm-1, nm-1);   // m x m
    // P is symmetric so P_mn = P_nm.t() — no separate copy needed

    // Split vecPos for quadratic-form computation
    const arma::vec vx = vecPos.subvec(0, n-1);
    const arma::vec vy = vecPos.subvec(n, nm-1);
    const double norm2 = arma::dot(vecPos, vecPos);
    // Note: vecNeg = [vx; -vy] has the same norm, so
    //   UNeg[l] = 0.5*(norm2 - 2*<vx, C_l*vy>) = norm2 - UPos[l]

    // Pre-flatten derivative matrices (column-major) for fast element-wise
    // trace: tr(A*B) = dot(vec(A), vec(B)) when B is symmetric
    arma::mat dX_flat(n * n, 3);
    arma::mat dY_flat(m * m, 3);
    for (arma::uword j = 0; j < 3; j++) {
        dX_flat.col(j) = arma::vectorise(derivX.slice(j));
        dY_flat.col(j) = arma::vectorise(derivY.slice(j));
    }

    arma::vec Itautau(L);
    arma::mat Itautheta(L, 6);
    arma::vec UPos(L), UNeg(L);

    for (arma::uword l = 0; l < L; l++) {
        // Non-copying view of slice l
        const arma::mat C(const_cast<double*>(crossBlocks.slice_memptr(l)),
                          n, m, false, true);

        // ---- PA_l = P * Sigma_alt_l = P + P * D_l ----
        // D_l = [[0, C]; [C^T, 0]], so P*D_l adds rank-2nm corrections:
        const arma::mat PA_nn = P_nn + P_nm * C.t();      // n x n
        const arma::mat PA_nm = P_nm + P_nn * C;          // n x m
        const arma::mat PA_mn = P_nm.t() + P_mm * C.t(); // m x n
        const arma::mat PA_mm = P_mm + P_nm.t() * C;     // m x m

        // Itautau[l] = 0.5 * ||PA_l||_F^2  (sum across all four blocks)
        Itautau(l) = 0.5 * (arma::accu(arma::square(PA_nn)) +
                             arma::accu(arma::square(PA_nm)) +
                             arma::accu(arma::square(PA_mn)) +
                             arma::accu(arma::square(PA_mm)));

        // ---- PA2 diagonal blocks: only needed for trace with derivX/Y ----
        // PA2_nn = (P * PA_l)[n-block] = P_nn*PA_nn + P_nm*PA_mn
        // PA2_mm = (P * PA_l)[m-block] = P_nm^T*PA_nm + P_mm*PA_mm
        const arma::mat PA2_nn = P_nn * PA_nn + P_nm * PA_mn;
        const arma::mat PA2_mm = P_nm.t() * PA_nm + P_mm * PA_mm;

        // Itautheta[l, j] = 0.5 * tr(PA2_nn * derivX_j) = 0.5 * dot(vec(PA2_nn), vec(derivX_j))
        const arma::vec pa2nn_v = arma::vectorise(PA2_nn);
        const arma::vec pa2mm_v = arma::vectorise(PA2_mm);
        for (arma::uword j = 0; j < 3; j++) {
            Itautheta(l, j)   = 0.5 * arma::dot(pa2nn_v, dX_flat.col(j));
            Itautheta(l, 3+j) = 0.5 * arma::dot(pa2mm_v, dY_flat.col(j));
        }

        // Score quadratic forms: U_l = 0.5 * (norm2 +/- 2 * <vx, C_l * vy>)
        const double cross = arma::dot(vx, C * vy);
        UPos(l) = 0.5 * (norm2 + 2.0 * cross);
        UNeg(l) = 0.5 * (norm2 - 2.0 * cross);
    }

    return Rcpp::List::create(
        Rcpp::Named("Itautau")   = Itautau,
        Rcpp::Named("Itautheta") = Itautheta,
        Rcpp::Named("UPos")      = UPos,
        Rcpp::Named("UNeg")      = UNeg
    );
}
