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

    // Batch: SigmaX * W_all in one DGEMM call (n x m*nWs)
    // Access the cube's contiguous memory as a flat matrix (no copy)
    const arma::mat W_mat(const_cast<double*>(W.memptr()), n, m * nWs, false, true);
    arma::mat SigW = SigmaX * W_mat;  // n x (m * nWs)

    arma::mat sigXws(mm2, nWs);
    // Use Rcpp::NumericVector so traces comes back as a plain R vector, not a matrix
    Rcpp::NumericVector traces(nWs);

    for (int wi = 0; wi < nWs; wi++) {
        // Views into pre-allocated memory (no copy)
        const arma::mat Wi (const_cast<double*>(W.slice_memptr(wi)),   n, m, false, true);
        const arma::mat SWi(SigW.colptr(wi * m),                       n, m, false, true);

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
