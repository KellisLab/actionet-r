// Rcpp interface for `tools` module
// Organized by module header in th order imported.
// [[Rcpp::interfaces(r, cpp)]]
#include "actionet_r_config.h"

// autocorrelation =====================================================================================================

// TODO: Unused. Remove?
// [[Rcpp::export]]
Rcpp::List
    autocorrelation_Moran_parametric(arma::sp_mat G, arma::mat scores, int normalization_method = 4,
                                     int thread_no = 0) {
    arma::field<arma::vec> out = actionet::autocorrelation_Moran_parametric(G, scores, normalization_method, thread_no);

    Rcpp::List res;
    res["stat"] = out[0];
    res["zscore"] = out[1];
    res["mu"] = out[2];
    res["sigma"] = out[3];

    return (res);
}

// [[Rcpp::export]]
Rcpp::List autocorrelation_Moran(arma::sp_mat G, arma::mat scores, int normalization_method = 1, int perm_no = 30,
                                 int thread_no = 0) {
    arma::field<arma::vec> out = actionet::autocorrelation_Moran(G, scores, normalization_method, perm_no, thread_no);

    Rcpp::List res;
    res["Moran_I"] = out[0];
    res["zscore"] = out[1];
    res["mu"] = out[2];
    res["sigma"] = out[3];

    return (res);
}

// [[Rcpp::export]]
Rcpp::List autocorrelation_Geary(arma::sp_mat G, arma::mat scores, int normalization_method = 1, int perm_no = 30,
                                 int thread_no = 0) {
    arma::field<arma::vec> out = actionet::autocorrelation_Geary(G, scores, normalization_method, perm_no, thread_no);

    Rcpp::List res;
    res["Geary_C"] = out[0];
    res["zscore"] = out[1];
    res["mu"] = out[2];
    res["sigma"] = out[3];

    return (res);
}

// enrichment ==========================================================================================================

// [[Rcpp::export]]
arma::mat assess_label_enrichment(arma::sp_mat& G, arma::mat& M, int thread_no = 0) {
    arma::mat logPvals = actionet::assess_label_enrichment(G, M, thread_no);

    return (logPvals);
}

//' Computes feature enrichment wrt a given annotation
//'
//' @param scores Specificity scores of features
//' @param associations Binary matrix of annotations
//' @param L Length of the top-ranked scores to scan
//'
//' @return Matrix of log-pvalues
//'
//' @examples
//' data("gProfilerDB_human")
//' G = colNets(ace)$ACTIONet
//' associations = gProfilerDB_human$SYMBOL$REAC
//' common.genes = intersect(rownames(ace), rownames(associations))
//' specificity_scores = rowFactors(ace)[["H_merged_upper_significance"]]
//' logPvals = compute_feature_specificity(
//' specificity_scores[common.genes, ], annotations[common.genes, ]
//' )
//' rownames(logPvals) = colnames(specificity_scores)
//' colnames(logPvals) = colnames(annotations)
// [[Rcpp::export]]
Rcpp::List assess_enrichment(arma::mat& scores, arma::sp_mat& associations, int thread_no = 0) {
    arma::field<arma::mat> res = actionet::assess_enrichment(scores, associations, thread_no);

    Rcpp::List out_list;
    out_list["logPvals"] = res(0);
    out_list["thresholds"] = res(1);

    return (out_list);
}

// matrix_misc =========================================================================================================

//' Aggregate matrix within groups
//'
//' @param S matrix of type "dMatrix"
//' @param sample_assignments Vector of column groupings. Group labels must be continuous integers or coercible to such.
//'
//' @return S matrix with columns of values aggregated within each group of sample_assignments
//'
// [[Rcpp::export]]
arma::mat compute_grouped_rowsums(arma::sp_mat& S, arma::vec sample_assignments) {
    arma::mat pb = actionet::compute_grouped_rowsums(S, sample_assignments);

    return pb;
}

// [[Rcpp::export]]
arma::mat compute_grouped_rowsums_full(arma::mat& S, arma::vec sample_assignments) {
    arma::mat pb = actionet::compute_grouped_rowsums(S, sample_assignments);

    return pb;
}

//' Average matrix within groups
//'
//' @param S matrix
//' @param sample_assignments Vector of column groupings. Group labels must be continuous integers or coercible to such.
//'
//' @return S matrix with columns of values average within each group of sample_assignments
//'
// [[Rcpp::export]]
arma::mat compute_grouped_rowmeans(arma::sp_mat& S, arma::vec sample_assignments) {
    arma::mat pb = actionet::compute_grouped_rowmeans(S, sample_assignments);

    return pb;
}

// [[Rcpp::export]]
arma::mat compute_grouped_rowmeans_full(arma::mat& S, arma::vec sample_assignments) {
    arma::mat pb = actionet::compute_grouped_rowmeans(S, sample_assignments);

    return pb;
}

// [[Rcpp::export]]
arma::mat compute_grouped_rowvars(arma::sp_mat& S, arma::vec sample_assignments) {
    arma::mat pb = actionet::compute_grouped_rowvars(S, sample_assignments);

    return pb;
}

// [[Rcpp::export]]
arma::mat compute_grouped_rowvars_full(arma::mat& S, arma::vec sample_assignments) {
    arma::mat pb = actionet::compute_grouped_rowvars(S, sample_assignments);

    return pb;
}

// mwm =================================================================================================================

//' Computes the maximum-weight bipartite graph matching
//'
//' @param G Adjacency matrix of the input graph
//'
//' @return G_matched An adjacency matrix with a maximum of one nonzero entry on
//' rows/columns
//'
//' @examples
//' G_matched = MWM_hungarian(G)
// [[Rcpp::export]]
arma::mat MWM_hungarian(arma::mat& G) {
    arma::mat G_matched = actionet::MWM_hungarian(G);

    return G_matched;
}

// [[Rcpp::export]]
arma::umat MWM_rank1(arma::vec u, arma::vec v, double u_threshold = 0, double v_threshold = 0) {
    arma::umat pairs = actionet::MWM_rank1(u, v, u_threshold, v_threshold);

    pairs = pairs + 1;

    return (pairs);
}

// normalization =======================================================================================================

// TODO: Update and remove. Single reference.
// [[Rcpp::export]]
arma::mat normalize_mat(arma::mat& X, int p = 0, int dim = 0) {
    arma::mat X_norm = actionet::normalize_matrix(X, p, dim);

    return (X_norm);
}

// TODO: Update and remove. Single reference.
// [[Rcpp::export]]
arma::sp_mat normalize_spmat(arma::sp_mat& X, int p = 0, int dim = 0) {
    arma::sp_mat X_norm = actionet::normalize_matrix(X, p, dim);

    return (X_norm);
}

// xicor ===============================================================================================================

// [[Rcpp::export]]
arma::vec xicor(arma::vec xvec, arma::vec yvec, bool compute_pval = true, int seed = 0) {
    arma::vec res = actionet::xicor(xvec, yvec, compute_pval, seed);

    return (res);
}

// [[Rcpp::export]]
Rcpp::List XICOR(arma::mat& X, arma::mat& Y, bool compute_pval = true, int seed = 0, int thread_no = 0) {
    arma::field<arma::mat> out = actionet::XICOR(X, Y, compute_pval, seed, thread_no);

    Rcpp::List res;
    res["XI"] = out(0);
    res["Z"] = out(1);

    return (res);
}
