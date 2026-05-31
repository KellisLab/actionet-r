// Rcpp interface for `visualization` module
// Organized by module header in th order imported.
#include "actionet_r_config.h"

// generate_layout =====================================================================================================

// [[Rcpp::export]]
arma::mat C_layoutNetwork(arma::sp_mat& G, arma::mat& initial_coordinates, std::string method = "umap",
                          unsigned int n_components = 2, float spread = 1, float min_dist = 1,
                          unsigned int n_epochs = 0,
                          float learning_rate = 1, float repulsion_strength = 1, float negative_sample_rate = 3,
                          bool approx_pow = true, bool pcg_rand = true, std::string rng_type = "",
                          bool batch = true, unsigned int grain_size = 1,
                          Rcpp::NumericVector ai = R_NilValue, Rcpp::NumericVector aj = R_NilValue,
                          int seed = 0, int thread_no = 0, bool verbose = true, float a = 0, float b = 0,
                          std::string opt_method = "adam", float alpha = -1, float beta1 = 0.5,
                          float beta2 = 0.9, float eps = 1e-7) {
    const std::size_t requested_threads = thread_no > 0 ? static_cast<std::size_t>(thread_no) : 0;
    OptimizerArgs opt_args(opt_method, alpha == -1 ? learning_rate : alpha, beta1, beta2, eps);
    UwotArgs uwot_args(
        method, n_components, spread, min_dist, n_epochs, learning_rate,
        repulsion_strength, negative_sample_rate, approx_pow, pcg_rand,
        batch, seed, requested_threads, grain_size, verbose, opt_args, rng_type
    );

    if (a != 0 || b != 0) {
        uwot_args.set_ab(a, b);
    }
    if (!ai.isNULL() && ai.size() > 0) {
        std::vector<float> ai_vec;
        ai_vec.reserve(ai.size());
        for (double value : ai) {
            ai_vec.push_back(static_cast<float>(value));
        }
        uwot_args.ai = std::move(ai_vec);
    }
    if (!aj.isNULL() && aj.size() > 0) {
        std::vector<float> aj_vec;
        aj_vec.reserve(aj.size());
        for (double value : aj) {
            aj_vec.push_back(static_cast<float>(value));
        }
        uwot_args.aj = std::move(aj_vec);
    }

    arma::mat coordinates = actionet::layoutNetwork(G, initial_coordinates, std::move(uwot_args));

    return coordinates;
}

// [[Rcpp::export]]
arma::mat C_computeNodeColors(const arma::mat& coordinates, int thread_no = 1) {
    arma::mat colors = actionet::computeNodeColors(coordinates, thread_no);
    return (colors);
}
