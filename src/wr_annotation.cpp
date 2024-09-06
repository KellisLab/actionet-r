// Rcpp interface for `annotation` module
// Organized by module header in th order imported.
// [[Rcpp::interfaces(r, cpp)]]
#include "actionet_r_config.h"

// marker_stats ========================================================================================================

// [[Rcpp::export]]
arma::mat compute_marker_aggregate_stats(arma::sp_mat& G, arma::sp_mat& S, arma::sp_mat& marker_mat,
                                         double alpha = 0.85, int max_it = 5, int thread_no = 0,
                                         bool ignore_baseline_expression = false) {
    arma::mat stats = actionet::compute_marker_aggregate_stats(G, S, marker_mat, alpha, max_it, thread_no,
                                                               ignore_baseline_expression);

    return (stats);
}

// [[Rcpp::export]]
Rcpp::List aggregate_genesets_vision(arma::sp_mat& G, arma::sp_mat& S, arma::mat& marker_mat,
                                     int norm_type = 2, double alpha = 0.85, int max_it = 5, double tol = 1E-8, int thread_no = 0) {
    arma::field<arma::mat> stats = actionet::aggregate_genesets_vision(G, S, marker_mat, norm_type,
                                                                       alpha, max_it, tol, thread_no);

    Rcpp::List res;

    res["stats_norm_smoothed"] = stats[0];
    res["stats_norm"] = stats[1];
    res["stats"] = stats[2];

    return (res);
}

// specificity =========================================================================================================

//' Compute feature specificity (from archetype footprints)
//'
//' @param S Input matrix (sparseMatrix)
//' @param H A soft membership matrix - Typically H_merged from the merge_archetypes() function.
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//' prune.out = collect_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = merge_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//' S.norm = renormalize_input_matrix(S, cell.clusters)
//' logPvals.list = compute_archetype_feature_specificity(S.norm, unification.out$H_merged)
//' specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
Rcpp::List compute_archetype_feature_specificity(arma::sp_mat& S, arma::mat& H, int thread_no = 0) {
    arma::field<arma::mat> res = actionet::compute_feature_specificity(S, H, thread_no);

    Rcpp::List out_list;
    out_list["archetypes"] = res(0);
    out_list["upper_significance"] = res(1);
    out_list["lower_significance"] = res(2);

    return (out_list);
}

// [[Rcpp::export]]
Rcpp::List compute_archetype_feature_specificity_full(arma::mat& S, arma::mat& H, int thread_no = 0) {
    arma::field<arma::mat> res = actionet::compute_feature_specificity(S, H, thread_no);

    Rcpp::List out_list;
    out_list["archetypes"] = res(0);
    out_list["upper_significance"] = res(1);
    out_list["lower_significance"] = res(2);

    return (out_list);
}

//' Compute feature specificity (from cluster assignments)
//'
//' @param S Input matrix ("sparseMatrix")
//' @param sample_assignments Vector of cluster assignments
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//' prune.out = collect_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = merge_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//' S.norm = renormalize_input_matrix(S, cell.clusters)
//' logPvals.list = compute_cluster_feature_specificity(S.norm, cell.clusters)
//' specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
Rcpp::List compute_cluster_feature_specificity(arma::sp_mat& S, arma::uvec sample_assignments, int thread_no = 0) {
    arma::field<arma::mat> res = actionet::compute_feature_specificity(S, sample_assignments, thread_no);

    Rcpp::List out_list;
    out_list["average_profile"] = res(0);
    out_list["upper_significance"] = res(1);
    out_list["lower_significance"] = res(2);

    return (out_list);
}

// [[Rcpp::export]]
Rcpp::List compute_cluster_feature_specificity_full(arma::mat& S, arma::uvec sample_assignments, int thread_no = 0) {
    arma::field<arma::mat> res = actionet::compute_feature_specificity(S, sample_assignments, thread_no);

    Rcpp::List out_list;
    out_list["average_profile"] = res(0);
    out_list["upper_significance"] = res(1);
    out_list["lower_significance"] = res(2);

    return (out_list);
}
