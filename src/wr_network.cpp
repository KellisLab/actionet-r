// Rcpp interface for `network` module
// Organized by module header in th order imported.
// [[Rcpp::interfaces(r, cpp)]]
#include "actionet_r_config.h"

// build_network =======================================================================================================

//' Builds an interaction network from the multi-level archetypal decompositions
//'
//' @param H_stacked Output of the collect_archetypes() function.
//' @param density Overall density of constructed graph. The higher the density,
//' the more edges are retained (default = 1.0).
//' @param thread_no Number of parallel threads (default = 0).
//' @param mutual_edges_only Symmetrization strategy for nearest-neighbor edges.
//' If it is true, only mutual nearest-neighbors are returned (default=TRUE).
//'
//' @return G Adjacency matrix of the ACTIONet graph.
//'
//' @examples
//' prune.out = collect_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
// [[Rcpp::export]]
arma::sp_mat buildNetwork(arma::mat H, std::string algorithm = "k*nn", std::string distance_metric = "jsd",
                          double density = 1.0, int thread_no = 0, bool mutual_edges_only = true, int k = 10) {
    // TODO: Add as options in interface
    double M = 16, ef_construction = 200, ef = 50;
    arma::sp_mat G = actionet::buildNetwork(H, algorithm, distance_metric, density, thread_no, M,
                                            ef_construction, ef, mutual_edges_only, k);

    return G;
}

// label_propagation ===================================================================================================

// [[Rcpp::export]]
arma::vec run_LPA(arma::sp_mat& G, arma::vec labels, double lambda = 1, int iters = 3,
                  double sig_threshold = 3, Rcpp::Nullable<Rcpp::IntegerVector> fixed_labels_ = R_NilValue,
                  int thread_no = 0) {
    arma::uvec fixed_labels_vec;
    if (fixed_labels_.isNotNull()) {
        Rcpp::NumericVector fixed_labels(fixed_labels_);
        fixed_labels_vec.set_size(fixed_labels.size());
        for (int i = 0; i < fixed_labels.size(); i++) {
            fixed_labels_vec(i) = fixed_labels(i) - 1;
        }
    }

    arma::vec new_labels =
        actionet::LPA(G, std::move(labels), lambda, iters, sig_threshold, fixed_labels_vec, thread_no);

    return (new_labels);
}

// network_diffusion ===================================================================================================

//' Computes network diffusion over a given network, starting with an arbitrarty
//' set of initial scores
//'
//' @param G Input graph
//' @param X0 Matrix of initial values per diffusion (ncol(G) == nrow(G) == ncol(X0))
//' @param thread_no Number of parallel threads (default=0)
//' @param alpha Random-walk depth ( between [0, 1] )
//' @param max_it PageRank iterations
//'
//' @return Matrix of diffusion scores
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' gene.expression = Matrix::t(logcounts(ace))[c("CD19", "CD14", "CD16"), ]
//' smoothed.expression = compute_network_diffusion(G, gene.expression)
// [[Rcpp::export]]
arma::mat compute_network_diffusion_fast(arma::sp_mat& G, arma::sp_mat& X0, double alpha = 0.85, int max_it = 3,
                                         int thread_no = 0) {
    arma::mat Diff = actionet::compute_network_diffusion_fast(G, X0, alpha, max_it, thread_no);

    return (Diff);
}

//' Computes network diffusion over a given network, starting with an arbitrarty
//' set of initial scores
//'
//' @param G Input graph
//' @param X0 Matrix of initial values per diffusion (ncol(G) == nrow(G) == ncol(X0))
//' @param thread_no Number of parallel threads (default=0)
//' @param alpha Random-walk depth ( between [0, 1] ) ' @param max_it PageRank iterations
//'
//' @return Matrix of diffusion scores
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' gene.expression = Matrix::t(logcounts(ace))[c("CD19", "CD14", "CD16"), ]
//' smoothed.expression = compute_network_diffusion_approx(G, gene.expression)
// [[Rcpp::export]]
arma::mat compute_network_diffusion_approx(arma::sp_mat& G, arma::mat& X0, int norm_type = 0, double alpha = 0.85,
                                           int max_it = 5, double tol = 1e-8, int thread_no = 0) {
    if (G.n_rows != X0.n_rows) {
        stderr_printf("Dimension mismatch: G (%dx%d) and X0 (%dx%d)\n", G.n_rows, G.n_cols, X0.n_rows, X0.n_cols);
        return (arma::mat());
    }

    arma::mat X = actionet::compute_network_diffusion_approx(G, X0, norm_type, alpha, max_it, tol, thread_no);

    return (X);
}

// network_measures ====================================================================================================

//' Compute coreness of graph vertices
//'
//' @param G Input graph
//'
//' @return cn core-number of each graph node
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' cn = compute_core_number(G)
// [[Rcpp::export]]
arma::uvec compute_core_number(arma::sp_mat& G) {
    arma::uvec core_num = actionet::compute_core_number(G);

    return (core_num);
}

//' Compute coreness of subgraph vertices induced by each archetype
//'
//' @param G Input graph
//' @param sample_assignments Archetype discretization (output of merge_archetypes())
//'
//' @return cn core-number of each graph node
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' assignments = ace$archetype.assignment
//' connectivity = compute_core_number(G, assignments)
// [[Rcpp::export]]
arma::vec compute_archetype_core_centrality(arma::sp_mat& G, arma::uvec sample_assignments) {
    arma::vec conn = actionet::compute_archetype_core_centrality(G, sample_assignments);

    return (conn);
}
