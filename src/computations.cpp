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
    const arma::cube& W,
    bool findSigXws
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
        if(findSigXws){
            arma::uword k = 0;
            for (int j = 0; j < m - 1; j++) {
                for (int i = j + 1; i < m; i++) {
                    sigXws(k++, wi) = tmp(i, j);
                }
            }
        }
        // Skip this step for point patterns
        traces[wi] = arma::trace(tmp);
    }

    return Rcpp::List::create(
        Rcpp::Named("sigXws") = sigXws,
        Rcpp::Named("traces") = traces
    );
}

// Helper: column-major pairwise Euclidean distances matching R's dist()
// i.e. for pairs (i,j) with j < i, iterating j from 0..n-2, i from j+1..n-1
static arma::vec pairwiseDist2D(const arma::mat& C) {
    arma::uword n   = C.n_rows;
    arma::uword nn2 = n * (n - 1) / 2;
    arma::vec   d(nn2);
    arma::uword idx = 0;
    for (arma::uword j = 0; j < n - 1; j++) {
        for (arma::uword i = j + 1; i < n; i++) {
            double dx = C(i, 0) - C(j, 0);
            double dy = C(i, 1) - C(j, 1);
            d[idx++] = std::sqrt(dx * dx + dy * dy);
        }
    }
    return d;
}

// Helper: first nn2 pairwise distances from the full m-point lower triangle
// in column-major order, stopping early once nn2 pairs are collected.
// Replicates dist(Ey)[1:nn2] in R when Ey has m rows — O(nn2) work, not O(mm2).
static arma::vec firstNDistsPPP(const arma::mat& Ey, arma::uword nn2) {
    arma::uword m = Ey.n_rows;
    arma::vec   d(nn2);
    arma::uword cnt = 0;
    for (arma::uword j = 0; j < m && cnt < nn2; j++) {
        for (arma::uword i = j + 1; i < m && cnt < nn2; i++) {
            double dx = Ey(i, 0) - Ey(j, 0);
            double dy = Ey(i, 1) - Ey(j, 1);
            d[cnt++] = std::sqrt(dx * dx + dy * dy);
        }
    }
    return d;
}

//' Variance traces for \code{MoransISinglePPP}: compute \eqn{tr(W^T \Sigma_{\!X} W)} for each Y feature
//'
//' In the PPP setting there is no fitted covariance model for the X modality, so
//' \eqn{\Sigma_{\!X}} is approximated by the \eqn{n \times n} covariance matrix
//' among the \emph{first n} Y spots (matching the existing \code{computeSigXws}
//' behaviour when called with \eqn{m(m-1)/2} Y-covariance values and an
//' \eqn{n \times m} weight matrix).  Distances are computed here from \code{Ey}
//' without materialising a distance vector in R, and the identity
//' \eqn{tr(W^T \Sigma_{\!X} W) = tr(\Sigma_{\!X} \cdot W W^T)} is exploited so
//' that, per feature, only \eqn{n \times n} arithmetic is needed instead of the
//' \eqn{n \times m} and \eqn{m \times m} intermediates in the original code.
//'
//' @param W      \eqn{n \times m} weight matrix (single slice, already normalised)
//' @param Ey     \eqn{m \times 2} coordinate matrix for the second modality
//' @param vgParY \eqn{k \times 3} variogram parameters for Y features:
//'   columns \code{[psill, range, isExp]}
//' @return Length-\eqn{k} vector of raw variance values (before division by
//'   \code{prodFac}); one entry per Y feature.
//' @keywords internal
// [[Rcpp::export]]
arma::vec computeTracePPP_cpp(
    const arma::mat& W,
    const arma::mat& Ey,
    const arma::mat& vgParY
) {
    arma::uword n   = W.n_rows;
    arma::uword k   = vgParY.n_rows;
    arma::uword nn2 = n * (n - 1) / 2;

    const arma::vec distFirst = firstNDistsPPP(Ey, nn2);

    arma::vec traces(k);
    arma::mat SigX(n, n);

    for (arma::uword fj = 0; fj < k; fj++) {
        arma::vec vgY = evalVariogramCpp(distFirst,
                                          vgParY(fj, 0),
                                          vgParY(fj, 1),
                                          vgParY(fj, 2) > 0.5);

        // Build SigmaX (n x n): diagonal = 1, off-diagonal from the first
        // nn2 variogram values in column-major lower-triangle order
        SigX.eye();
        {
            arma::uword idx = 0;
            for (arma::uword j2 = 0; j2 < n - 1; j2++) {
                for (arma::uword i2 = j2 + 1; i2 < n; i2++) {
                    SigX(i2, j2) = vgY[idx];
                    SigX(j2, i2) = vgY[idx];
                    idx++;
                }
            }
        }

        // tr(W^T SigX W) — same computation as one slice of computeSigXws
        arma::mat SWi = SigX * W;                  // n x m
        traces(fj) = arma::trace(W.t() * SWi);    // scalar via m x m
    }

    return traces;
}

//' Build the Gaussian weight matrix and compute Ixy + variance traces for MoransISinglePPP
//'
//' Combines weight-matrix construction, Ixy calculation, and variance-trace
//' computation in a single C++ call so that the \eqn{n \times m} weight matrix
//' \eqn{W} is never materialised as an R object.  The weight matrix is
//' \eqn{W_{ij} \propto \exp(-\|C_x^{(i)} - E_y^{(j)}\|^2 / \eta)}, normalised
//' to sum to one.
//'
//' @param Cx          \eqn{n \times 2} X-point coordinates for the current cell type
//' @param Ey          \eqn{m \times 2} Y-spot coordinates
//' @param eta         Gaussian bandwidth \eqn{\eta}
//' @param Y           \eqn{m \times k} scaled Y feature matrix (columns = \code{featuresY})
//' @param vgParY      \eqn{k \times 3} variogram parameters \code{[psill, range, isExp]}
//'   (ignored when \code{findVariances = FALSE})
//' @param sqrtProdFac \eqn{\sqrt{(n-1)(m-1)}} normalisation factor
//' @param findVariances logical; whether to compute variances
//' @return A list with
//'   \describe{
//'     \item{isZero}{logical; \code{TRUE} if \eqn{W} sums to zero (all weights underflow)}
//'     \item{Ixys}{length-\eqn{k} vector of Ixy statistics (zero when \code{isZero})}
//'     \item{traces}{length-\eqn{k} vector of \eqn{tr(W^T \Sigma_X W)} (zero when
//'       \code{!findVariances} or \code{isZero})}
//'     \item{trWtW}{\eqn{tr(W^T W) = \sum W_{ij}^2}, used as the independence fallback}
//'   }
//' @keywords internal
//' @note This function is highly optmised to keep memory usage low, at an elevated computation cost
// [[Rcpp::export]]
Rcpp::List computeIxyAndTracePPP_cpp(
        const arma::mat& Cx,
        const arma::mat& Ey,
        double           eta,
        const arma::mat& Y,
        const arma::mat& vgParY,
        double           sqrtProdFac,
        bool             findVariances
) {
    arma::uword n = Cx.n_rows;
    arma::uword m = Ey.n_rows;
    arma::uword k = Y.n_cols;

    // Instead of building the full (n x m) W, compute one column of W at a
    // time, accumulate what's needed, and discard the column immediately.
    arma::vec wColSums(m, arma::fill::zeros);   // raw (unnormalised) column sums of W
    arma::mat WWT(n, n, arma::fill::zeros);     // raw (unnormalised) W * W^T, built incrementally
    double Wsum = 0.0;

    arma::vec wcol(n);
    for (arma::uword j = 0; j < m; j++) {
        for (arma::uword i = 0; i < n; i++) {
            double dx = Cx(i, 0) - Ey(j, 0);
            double dy = Cx(i, 1) - Ey(j, 1);
            wcol(i) = std::exp(-(dx * dx + dy * dy) / eta);
        }
        double colSum = arma::sum(wcol);
        wColSums(j) = colSum;
        Wsum += colSum;
        if (findVariances) {
            WWT += wcol * wcol.t();   // outer product contributes to W * W^T
        }
        // wcol is reused/overwritten next iteration -- never stored as a full matrix
    }

    if (Wsum == 0.0) {
        return Rcpp::List::create(
            Rcpp::Named("isZero") = true,
            Rcpp::Named("Ixys")   = arma::zeros<arma::vec>(k),
            Rcpp::Named("traces") = arma::zeros<arma::vec>(k),
            Rcpp::Named("trWtW")  = 0.0
        );
    }

    // Normalise at the end.
    wColSums /= Wsum;
    const arma::vec Ixys = (Y.t() * wColSums) / sqrtProdFac;  // k-vector

    if (!findVariances) {
        return Rcpp::List::create(
            Rcpp::Named("isZero") = false,
            Rcpp::Named("Ixys")   = Ixys,
            Rcpp::Named("traces") = arma::zeros<arma::vec>(k),
            Rcpp::Named("trWtW")  = 0.0
        );
    }

    // WWT was accumulated from the raw (unnormalised) W; normalise now.
    WWT /= (Wsum * Wsum);
    const double trWtW = arma::trace(WWT);

    const arma::uword nn2 = n * (n - 1) / 2;
    const arma::vec distFirst = firstNDistsPPP(Ey, nn2);
    arma::mat SigX(n, n);
    arma::vec traces(k);
    for (arma::uword fj = 0; fj < k; fj++) {
        arma::vec vgY = evalVariogramCpp(distFirst,
                                         vgParY(fj, 0), vgParY(fj, 1),
                                         vgParY(fj, 2) > 0.5);
        SigX.eye();
        {
            arma::uword idx = 0;
            for (arma::uword j2 = 0; j2 < n - 1; j2++) {
                for (arma::uword i2 = j2 + 1; i2 < n; i2++) {
                    SigX(i2, j2) = vgY[idx];
                    SigX(j2, i2) = vgY[idx];
                    idx++;
                }
            }
        }
        traces(fj) = arma::trace(SigX * WWT);
    }
    return Rcpp::List::create(
        Rcpp::Named("isZero") = false,
        Rcpp::Named("Ixys")   = Ixys,
        Rcpp::Named("traces") = traces,
        Rcpp::Named("trWtW")  = trWtW
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
