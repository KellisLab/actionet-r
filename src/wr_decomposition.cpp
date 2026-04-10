// Rcpp interface for `decomposition` module
// Organized by module header in th order imported.
//
// NOTE (Plan 02): After the C++ core contract flip, the expected orientations
// for this reference wrapper copy are:
//   S          — cells x genes  (obs x var, AnnData-native)
//   S_r (in)   — cells x k
//   old_S_r    — cells x k  (public reduction contract; pass directly)
//   design     — cells x covariates  (unchanged)
// The actual actionet-r wrappers will be updated in Plan 03.
#include "actionet_r_config.h"

// orthogonalization ==========================================================================================================
// [[Rcpp::export]]
Rcpp::List C_orthogonalizeBatchEffect(arma::sp_mat& S, arma::mat& old_S_r, arma::mat& old_U, arma::mat& old_A,
                                      arma::mat& old_B, arma::vec& old_sigma, arma::mat& design) {
    arma::field<arma::mat> reduction_results(5);
    reduction_results(0) = old_S_r;
    reduction_results(1) = old_sigma;
    reduction_results(2) = old_U;
    reduction_results(3) = old_A;
    reduction_results(4) = old_B;

    arma::field<arma::mat> orthogonalized_reduction = actionet::orthogonalizeBatchEffect(S, reduction_results, design);

    Rcpp::List res;
    res["S_r"] = orthogonalized_reduction(0);

    arma::vec sigma = orthogonalized_reduction(1).col(0);
    res["sigma"] = sigma;
    res["U"] = orthogonalized_reduction(2);

    res["A"] = orthogonalized_reduction(3);
    res["B"] = orthogonalized_reduction(4);

    return res;
}

//[[Rcpp::export]]
Rcpp::List C_orthogonalizeBatchEffect_full(arma::mat& S, arma::mat& old_S_r, arma::mat& old_U, arma::mat& old_A,
                                           arma::mat& old_B, arma::vec& old_sigma, arma::mat& design) {
    arma::field<arma::mat> reduction_results(5);
    reduction_results(0) = old_S_r;
    reduction_results(1) = old_sigma;
    reduction_results(2) = old_U;
    reduction_results(3) = old_A;
    reduction_results(4) = old_B;

    arma::field<arma::mat> orthogonalized_reduction = actionet::orthogonalizeBatchEffect(S, reduction_results, design);

    Rcpp::List res;
    res["S_r"] = orthogonalized_reduction(0);

    arma::vec sigma = orthogonalized_reduction(1).col(0);
    res["sigma"] = sigma;
    res["U"] = orthogonalized_reduction(2);
    res["A"] = orthogonalized_reduction(3);
    res["B"] = orthogonalized_reduction(4);

    return res;
}

// [[Rcpp::export]]
Rcpp::List C_orthogonalizeBasal(arma::sp_mat& S, arma::mat& old_S_r, arma::mat& old_U, arma::mat& old_A,
                                arma::mat& old_B, arma::vec& old_sigma, arma::mat& basal) {
    arma::field<arma::mat> reduction_results(5);
    reduction_results(0) = old_S_r;
    reduction_results(1) = old_sigma;
    reduction_results(2) = old_U;
    reduction_results(3) = old_A;
    reduction_results(4) = old_B;

    arma::field<arma::mat> orthogonalized_reduction = actionet::orthogonalizeBasal(S, reduction_results, basal);

    Rcpp::List res;
    res["S_r"] = orthogonalized_reduction(0);

    arma::vec sigma = orthogonalized_reduction(1).col(0);
    res["sigma"] = sigma;
    res["U"] = orthogonalized_reduction(2);

    res["A"] = orthogonalized_reduction(3);
    res["B"] = orthogonalized_reduction(4);

    return res;
}

//[[Rcpp::export]]
Rcpp::List C_orthogonalizeBasal_full(arma::mat& S, arma::mat& old_S_r, arma::mat& old_U, arma::mat& old_A,
                                     arma::mat& old_B, arma::vec& old_sigma, arma::mat& basal) {
    arma::field<arma::mat> reduction_results(5);
    reduction_results(0) = old_S_r;
    reduction_results(1) = old_sigma;
    reduction_results(2) = old_U;
    reduction_results(3) = old_A;
    reduction_results(4) = old_B;

    arma::field<arma::mat> orthogonalized_reduction = actionet::orthogonalizeBasal(S, reduction_results, basal);

    Rcpp::List res;
    res["S_r"] = orthogonalized_reduction(0);

    arma::vec sigma = orthogonalized_reduction(1).col(0);
    res["sigma"] = sigma;
    res["U"] = orthogonalized_reduction(2);
    res["A"] = orthogonalized_reduction(3);
    res["B"] = orthogonalized_reduction(4);

    return res;
}

// svd_main =========================================================================================================

//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm:
//' From: IRLBA R Package
//'
//' @param A Input matrix ("sparseMatrix")
//' @param k Dimension of SVD decomposition
//' @param max_it Number of iterations (default=5)
//' @param seed Random seed (default=0)
//' @param algorithm SVD algorithm to use:
//'   - 0 = IRLB (default, good for small/medium matrices)
//'   - 1 = Halko (randomized SVD)
//'   - 2 = Feng (another randomized method)
//'   - 3 = PRIMME (NOT AVAILABLE in R - use Python for large matrices)
//'
//' @return A named list with U, sigma, and V components
//'
//' @note PRIMME algorithm is not available in R builds due to R's 32-bit matrix limitations.
//'   For very large sparse matrices (>2^31 elements), use the Python package.
//'
//' @examples
//' A = randn(100, 20)
//' svd.out = runSVD(A, dim = 3)
//' U = svd.out$u
// [[Rcpp::export]]
Rcpp::List C_runSVDSparse(arma::sp_mat& A, int k = 30, int max_it = 0, int seed = 0, int algorithm = 0,
                          bool verbose = true) {
#ifdef LIBACTIONET_BUILD_R
    // Guard against PRIMME usage in R builds
    if (algorithm == 3) {
        Rcpp::stop("PRIMME algorithm (3) is not available in R builds. R is limited to 32-bit matrix indices.\nFor large matrices, use the Python package with PRIMME support.");
    }
#endif

    arma::field<arma::mat> SVD_out = actionet::runSVD(A, k, max_it, seed, algorithm, verbose);

    Rcpp::List res;
    res["u"] = SVD_out(0);
    res["d"] = SVD_out(1);
    res["v"] = SVD_out(2);

    return res;
}

// [[Rcpp::export]]
Rcpp::List C_runSVDDense(arma::mat& A, int k = 30, int max_it = 0, int seed = 0, int algorithm = 0,
                         bool verbose = true) {
#ifdef LIBACTIONET_BUILD_R
    // Guard against PRIMME usage in R builds
    if (algorithm == 3) {
        Rcpp::stop("PRIMME algorithm (3) is not available in R builds. R is limited to 32-bit matrix indices.\nFor large matrices, use the Python package with PRIMME support.");
    }
#endif

    arma::field<arma::mat> SVD_out = actionet::runSVD(A, k, max_it, seed, algorithm, verbose);

    Rcpp::List res;
    res["u"] = SVD_out(0);
    res["d"] = SVD_out(1);
    res["v"] = SVD_out(2);

    return res;
}

// [[Rcpp::export]]
Rcpp::List C_perturbedSVD(const arma::mat& u, const arma::vec& d, const arma::mat& v,
                          const arma::mat& A, const arma::mat& B) {
    arma::vec d_vec = (d.n_cols > 1) ? arma::vec(d.diag()) : d;

    arma::field<arma::mat> SVD_results(3);
    SVD_results(0) = u;
    SVD_results(1) = d_vec;
    SVD_results(2) = v;

    arma::field<arma::mat> perturbed_SVD = actionet::perturbedSVD(SVD_results, A, B);

    Rcpp::List res;
    res["u"] = perturbed_SVD(0);
    res["d"] = perturbed_SVD(1).col(0);
    res["v"] = perturbed_SVD(2);

    return res;
}
