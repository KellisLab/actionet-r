// Rcpp interface for `visualization` module
// Organized by module header in th order imported.
// [[Rcpp::interfaces(r, cpp)]]
#include "actionet_r_config.h"

// generate_layout =====================================================================================================

//' Performs stochastic force-directed layout on the input graph (ACTIONet)
//'
//' @param G Adjacency matrix of the ACTIONet graph
//' @param S_r Reduced kernel matrix (is used for reproducible initialization).
//' @param compactness_level A value between 0-100, indicating the compactness
//' of ACTIONet layout (default=50)
//' @param n_epochs Number of epochs for SGD algorithm (default=100).
//' @param thread_no Number of threads (default = 0).
//'
//' @return A named list \itemize{
//' \item coordinates 2D coordinates of vertices.
//' \item coordinates_3D 3D coordinates of vertices.
//' \item colors De novo color of nodes inferred from their 3D embedding.
//' }
//'
//' @examples
//'	G = buildNetwork(prune.out$H_stacked)
//'	vis.out = layoutNetwrok(G, S_r)
// [[Rcpp::export]]
Rcpp::List layoutNetwork(arma::sp_mat& G, arma::mat& initial_position, const std::string& method = "umap",
                         double min_dist = 1, double spread = 1, double gamma = 1.0, unsigned int n_epochs = 500,
                         double learning_rate = 1.0, int seed = 0, int thread_no = 0) {
    arma::field<arma::mat> res = actionet::layoutNetwork_xmap(G, initial_position, method, min_dist, spread, gamma,
                                                              n_epochs, learning_rate, seed, thread_no);

    Rcpp::List out_list;
    out_list["coordinates"] = res(0);
    out_list["coordinates_3D"] = res(1);
    out_list["colors"] = res(2);

    return out_list;
}
