// Rcpp interface for libactionet
// [[Rcpp::interfaces(r, cpp)]]

// Enable build configuration for R interface
#define LIBACTIONET_BUILD_R
// Header `libactionet.hpp` configures package and includes `RcppArmadillo.h`. It must precede `Rcpp.h`.
#include "libactionet.hpp"

// [[Rcpp::depends(RcppArmadillo)]]

// TODO: Necessary?
template <typename T>
Rcpp::NumericVector arma2vec(const T& x)
{
    return Rcpp::NumericVector(x.begin(), x.end());
}

//' Set the RNG Seed from within Rcpp
//'
//' Within Rcpp, one can set the R session seed without triggering
//' the CRAN rng modifier check.
//' @param seed A \code{double} that is the seed one wishes to use.
//' @return A set RNG scope.
//' @examples
//' set.seed(10)
//' x = rnorm(5,0,1)
//' set_seed(10)
//' y = rnorm(5,0,1)
//' all.equal(x,y, check.attributes = F)
// [[Rcpp::export]]
void set_seed(double seed)
{
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}

//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm:
//' From: IRLBA R Package
//'
//' @param A Input matrix ("sparseMatrix")
//' @param k Dimension of SVD decomposition
//' @param max_it Number of iterations (default=5)
//' @param seed Random seed (default=0)
//' @param algorithm SVD algorithm to use. Currently supported methods are blah blah blah
//'
//' @return A named list with U, sigma, and V components
//'
//' @examples
//' A = randn(100, 20)
//' svd.out = runSVD(A, dim = 3)
//' U = svd.out$u
// [[Rcpp::export]]
Rcpp::List runSVD(arma::sp_mat& A, int k = 30, int max_it = 0, int seed = 0, int algorithm = 0, int verbose = 1)
{
    arma::field<arma::mat> SVD_out = actionet::runSVD(A, k, max_it, seed, algorithm, verbose);

    Rcpp::List res;
    res["u"] = SVD_out(0);
    res["d"] = SVD_out(1);
    res["v"] = SVD_out(2);

    return res;
}

// [[Rcpp::export]]
Rcpp::List runSVD_full(arma::mat& A, int k = 30, int max_it = 0, int seed = 0, int algorithm = 0, int verbose = 1)
{
    arma::field<arma::mat> SVD_out = actionet::runSVD(A, k, max_it, seed, algorithm, verbose);

    Rcpp::List res;
    res["u"] = SVD_out(0);
    res["d"] = SVD_out(1);
    res["v"] = SVD_out(2);

    return res;
}

// [[Rcpp::export]]
Rcpp::List perturbedSVD(arma::mat u, arma::vec d, arma::mat v, arma::mat A, arma::mat B)
{
    if (1 < d.n_cols)
        d = d.diag();

    arma::field<arma::mat> SVD_results(3);
    SVD_results(0) = u;
    SVD_results(1) = d;
    SVD_results(2) = v;

    arma::field<arma::mat> perturbed_SVD = actionet::perturbedSVD(SVD_results, A, B);

    Rcpp::List res;
    res["u"] = perturbed_SVD(0);
    res["d"] = perturbed_SVD(1).col(0);
    res["v"] = perturbed_SVD(2);

    return res;
}

// [[Rcpp::export]]
Rcpp::List orthogonalize_batch_effect(arma::sp_mat& S, arma::mat& old_S_r, arma::mat& old_V, arma::mat& old_A,
                                      arma::mat& old_B, arma::vec& old_sigma, arma::mat& design)
{
    arma::field<arma::mat> SVD_results(5);

    SVD_results(0) = old_V;
    SVD_results(1) = old_sigma;
    SVD_results(2) = old_S_r;
    for (int i = 0; i < old_sigma.n_elem; i++)
    {
        SVD_results(2).col(i) /= old_sigma(i);
    }
    SVD_results(3) = old_A;
    SVD_results(4) = old_B;

    arma::field<arma::mat> orthogonalized_reduction =
        actionet::orthogonalize_batch_effect(S, SVD_results, design);

    Rcpp::List res;
    res["V"] = orthogonalized_reduction(0);

    arma::vec sigma = orthogonalized_reduction(1).col(0);
    res["sigma"] = sigma;

    arma::mat V = orthogonalized_reduction(2);
    for (int i = 0; i < V.n_cols; i++)
    {
        V.col(i) *= sigma(i);
    }
    res["S_r"] = arma::trans(V);

    res["A"] = orthogonalized_reduction(3);
    res["B"] = orthogonalized_reduction(4);

    return res;
}

//[[Rcpp::export]]
Rcpp::List orthogonalize_batch_effect_full(arma::mat& S, arma::mat& old_S_r, arma::mat& old_V, arma::mat& old_A,
                                           arma::mat& old_B, arma::vec& old_sigma, arma::mat& design)
{
    arma::field<arma::mat> SVD_results(5);

    SVD_results(0) = old_V;
    SVD_results(1) = old_sigma;
    SVD_results(2) = old_S_r;
    for (int i = 0; i < old_sigma.n_elem; i++)
    {
        SVD_results(2).col(i) /= old_sigma(i);
    }
    SVD_results(3) = old_A;
    SVD_results(4) = old_B;

    arma::field<arma::mat> orthogonalized_reduction =
        actionet::orthogonalize_batch_effect(S, SVD_results, design);

    Rcpp::List res;
    res["V"] = orthogonalized_reduction(0);

    arma::vec sigma = orthogonalized_reduction(1).col(0);
    res["sigma"] = sigma;

    arma::mat V = orthogonalized_reduction(2);
    for (int i = 0; i < V.n_cols; i++)
    {
        V.col(i) *= sigma(i);
    }
    res["S_r"] = arma::trans(V);
    res["A"] = orthogonalized_reduction(3);
    res["B"] = orthogonalized_reduction(4);

    return res;
}

// [[Rcpp::export]]
Rcpp::List orthogonalize_basal(arma::sp_mat& S, arma::mat& old_S_r, arma::mat& old_V, arma::mat& old_A,
                               arma::mat& old_B, arma::vec& old_sigma, arma::mat& basal)
{
    arma::field<arma::mat> SVD_results(5);

    SVD_results(0) = old_V;
    SVD_results(1) = old_sigma;
    SVD_results(2) = old_S_r;
    for (int i = 0; i < old_sigma.n_elem; i++)
    {
        SVD_results(2).col(i) /= old_sigma(i);
    }
    SVD_results(3) = old_A;
    SVD_results(4) = old_B;

    arma::field<arma::mat> orthogonalized_reduction = actionet::orthogonalize_basal(S, SVD_results, basal);

    Rcpp::List res;
    res["V"] = orthogonalized_reduction(0);

    arma::vec sigma = orthogonalized_reduction(1).col(0);
    res["sigma"] = sigma;

    arma::mat V = orthogonalized_reduction(2);
    for (int i = 0; i < V.n_cols; i++)
    {
        V.col(i) *= sigma(i);
    }
    res["S_r"] = arma::trans(V);

    res["A"] = orthogonalized_reduction(3);
    res["B"] = orthogonalized_reduction(4);

    return res;
}

//[[Rcpp::export]]
Rcpp::List orthogonalize_basal_full(arma::mat& S, arma::mat& old_S_r, arma::mat& old_V, arma::mat& old_A,
                                    arma::mat& old_B, arma::vec& old_sigma, arma::mat& basal)
{
    arma::field<arma::mat> SVD_results(5);

    SVD_results(0) = old_V;
    SVD_results(1) = old_sigma;
    SVD_results(2) = old_S_r;
    for (int i = 0; i < old_sigma.n_elem; i++)
    {
        SVD_results(2).col(i) /= old_sigma(i);
    }
    SVD_results(3) = old_A;
    SVD_results(4) = old_B;

    arma::field<arma::mat> orthogonalized_reduction = actionet::orthogonalize_basal(S, SVD_results, basal);

    Rcpp::List res;
    res["V"] = orthogonalized_reduction(0);

    arma::vec sigma = orthogonalized_reduction(1).col(0);
    res["sigma"] = sigma;

    arma::mat V = orthogonalized_reduction(2);
    for (int i = 0; i < V.n_cols; i++)
    {
        V.col(i) *= sigma(i);
    }
    res["S_r"] = arma::trans(V);
    res["A"] = orthogonalized_reduction(3);
    res["B"] = orthogonalized_reduction(4);

    return res;
}

//' Computes reduced kernel matrix for a given (single-cell) profile
//'
//' @param S Input matrix ("sparseMatrix")
//' @param reduced_dim Dimension of the reduced kernel matrix (default=50)
//' @param iters Number of SVD iterations (default=5)
//' @param seed Random seed (default=0)
//' @param reduction_algorithm Kernel reduction algorithm. Currently only ACTION
//' method (1) is implemented (default=1)
//' @param SVD_algorithm SVD algorithm to use. Currently supported methods
//' are Halko (1) and Feng (2) (default=1)
//'
//' @return A named list with S_r, V, lambda, and exp_var. \itemize{
//' \item S_r: reduced kernel matrix of size reduced_dim x #samples.
//' \item V: Associated left singular-vectors (useful for reconstructing
//' discriminative scores for features, such as genes).
//' \item lambda, exp_var: Summary statistics of the sigular-values.
//' }
//'
//' @examples
//' S = logcounts(sce)
//' reduction.out = reduce(S, reduced_dim = 50)
//' S_r = reduction.out$S_r
// [[Rcpp::export]]
Rcpp::List reduce_kernel(arma::sp_mat& S, int reduced_dim = 50, int iter = 5, int seed = 0,
                         int SVD_algorithm = 0, bool prenormalize = false, int verbose = 1)
{
    arma::field<arma::mat> reduction =
        actionet::reduce_kernel(S, reduced_dim, iter, seed, SVD_algorithm, prenormalize, verbose);

    Rcpp::List res;
    res["V"] = reduction(0);

    arma::vec sigma = reduction(1).col(0);
    res["sigma"] = sigma;

    double epsilon = 0.01 / std::sqrt(reduction(2).n_rows);
    arma::mat V = arma::round(reduction(2) / epsilon) * epsilon;

    for (int i = 0; i < V.n_cols; i++)
    {
        arma::vec v = V.col(i) * sigma(i);
        V.col(i) = v;
    }
    V = arma::trans(V);
    res["S_r"] = V.eval();

    res["A"] = reduction(3);
    res["B"] = reduction(4);

    return res;
}

//' Computes reduced kernel matrix for a given (single-cell) profile
//'
//' @param S Input matrix ("matrix")
//' @param reduced_dim Dimension of the reduced kernel matrix (default=50)
//' @param iters Number of SVD iterations (default=5)
//' @param seed Random seed (default=0)
//' @param reduction_algorithm Kernel reduction algorithm. Currently only ACTION
//' method (1) is implemented (default=1)
//' @param SVD_algorithm SVD algorithm to use. Currently supported methods are
//' Halko (1) and Feng (2) (default=1)
//'
//' @return A named list with S_r, V, lambda, and exp_var. \itemize{
//' \item S_r: reduced kernel matrix of size reduced_dim x #samples.
//' \item V: Associated left singular-vectors (useful for reconstructing
//' discriminative scores for features, such as genes).
//' \item lambda, exp_var: Summary statistics of the sigular-values.
//' }
//'
//' @examples
//' S = logcounts(sce)
//' reduction.out = reduce(S, reduced_dim = 50)
//' S_r = reduction.out$S_r
// [[Rcpp::export]]
Rcpp::List reduce_kernel_full(arma::mat& S, int reduced_dim = 50, int iter = 5, int seed = 0, int SVD_algorithm = 0,
                              bool prenormalize = false, int verbose = 1)
{
    arma::field<arma::mat> reduction =
        actionet::reduce_kernel(S, reduced_dim, iter, seed, SVD_algorithm, prenormalize, verbose);

    Rcpp::List res;
    res["V"] = reduction(0);

    arma::vec sigma = reduction(1).col(0);
    res["sigma"] = sigma;

    double epsilon = 0.01 / std::sqrt(reduction(2).n_rows);
    arma::mat V = arma::round(reduction(2) / epsilon) * epsilon;

    for (int i = 0; i < V.n_cols; i++)
    {
        arma::vec v = V.col(i) * sigma(i);
        V.col(i) = v;
    }
    V = arma::trans(V);
    res["S_r"] = V.eval();

    res["A"] = reduction(3);
    res["B"] = reduction(4);

    return res;
}

//' Solves min_{X} (|| AX - B ||) s.t. simplex constraint
//'
//' @param A Input matrix
//' @param B Input matrix
//'
//' @return X Solution
//'
//' @examples
//' C = ACTION.out$C[[10]]
//' A = S_r %*% C
//' B = S_r
//' H = run_simplex_regression(A, B)
// [[Rcpp::export]]
arma::mat run_simplex_regression(arma::mat& A, arma::mat& B, bool computeXtX = false)
{
    arma::mat X = actionet::run_simplex_regression(A, B, computeXtX);

    return X;
}

//' Runs Successive Projection Algorithm (SPA) to solve separable NMF
//'
//' @param A Input matrix
//' @param k Number of columns to select
//'
//' @return A named list with entries 'selected_columns' and 'norms'
//' @examples
//' H = run_SPA(S_r, 10)
// [[Rcpp::export]]
Rcpp::List run_SPA(arma::mat& A, int k)
{
    actionet::SPA_results res = actionet::run_SPA(A, k);
    arma::uvec selected_columns = res.selected_columns;

    arma::vec cols(k);
    for (int i = 0; i < k; i++)
    {
        cols[i] = selected_columns[i] + 1;
    }

    Rcpp::List out;
    out["selected_columns"] = cols;
    out["norms"] = res.column_norms;

    return out;
}

//' Runs multi-level ACTION decomposition method
//'
//' @param S_r Reduced kernel matrix
//' @param k_min Minimum number of archetypes to consider (default=2)
//' @param k_max Maximum number of archetypes to consider, or "depth" of
//' decomposition (default=30)
//' @param thread_no Number of parallel threads (default = 0)
//' @param max_it,min_delta Convergence parameters for archetypal analysis
//'
//' @return A named list with entries 'C' and 'H', each a list for different
//' values of k ' @examples ' ACTION.out = run_ACTION(S_r, k_max = 10) ' H8 =
//' ACTION.out$H[[8]] ' cell.assignments = apply(H8, 2, which.max)
// [[Rcpp::export]]
Rcpp::List run_ACTION(arma::mat& S_r, int k_min = 2, int k_max = 30, int normalization = 1, int max_it = 100,
                      double min_delta = 1e-6, int thread_no = 0)
{
    actionet::ACTION_results trace =
        actionet::run_ACTION(S_r, k_min, k_max, normalization, max_it, min_delta, thread_no);

    Rcpp::List res;

    Rcpp::List C(k_max);
    for (int i = k_min; i <= k_max; i++)
    {
        arma::mat cur_C = trace.C[i];
        C[i - 1] = cur_C;
    }
    res["C"] = C;

    Rcpp::List H(k_max);
    for (int i = k_min; i <= k_max; i++)
    {
        arma::mat cur_H = trace.H[i];
        H[i - 1] = cur_H;
    }
    res["H"] = H;

    return res;
}

//' Runs basic archetypal analysis
//'
//' @param A Inpu matrix
//' @param W0 Starting archetypes
//' @param max_it,min_delta Convergence parameters for archetypal analysis
//'
//' @return A named list with entries 'C' and 'H', each a list for different
//' values of k ' @examples ' S_r = t(reducedDims(ace)$ACTION) ' SPA.out =
//' run_SPA(S_r, 10) ' W0 = S_r[, SPA.out$selected_columns] ' AA.out =
//' run_AA(S_r, W0) ' H = AA.out$H ' cell.assignments = apply(H, 2, which.max)
// [[Rcpp::export]]
Rcpp::List run_AA(arma::mat& A, arma::mat& W0, int max_it = 100, double min_delta = 1e-6)
{
    arma::field<arma::mat> res = actionet::run_AA(A, W0, max_it, min_delta);

    Rcpp::List out;
    out["C"] = res(0);
    out["H"] = res(1);
    out["W"] = A * res(0);

    return out;
}

//' Filters multi-level archetypes and concatenate filtered archetypes.
//' (Pre-ACTIONet archetype processing)
//'
//' @param C_trace,H_trace Output of ACTION
//' @param min_specificity_z_threshold Defines the stringency of pruning
//' nonspecific archetypes. The larger the value, the more archetypes will be
//' filtered out (default=-1)
//'
//' @return A named list: \itemize{
//' \item selected_archs: List of final archetypes that passed the
//' filtering/pruning step.
//' \item C_stacked,H_stacked: Horizontal/Vertical
//' concatenation of filtered C and H matrices, respectively.
//' }
//'
//' @examples
//' S = logcounts(sce)
//' reduction.out = reduce(S, reduced_dim = 50)
//' S_r = reduction.out$S_r
//' ACTION.out = run_ACTION(S_r, k_max = 10)
//' reconstruction.out = reconstruct_archetypes(S, ACTION.out$C, ACTION.out$H)
// [[Rcpp::export]]
Rcpp::List prune_archetypes(const Rcpp::List& C_trace, const Rcpp::List& H_trace,
                            double min_specificity_z_threshold = -3, int min_cells = 3)
{
    int n_list = H_trace.size();
    arma::field<arma::mat> C_trace_vec(n_list + 1);
    arma::field<arma::mat> H_trace_vec(n_list + 1);
    for (int i = 0; i < n_list; i++)
    {
        if (Rf_isNull(H_trace[i]))
        {
            continue;
        }
        C_trace_vec[i + 1] = (Rcpp::as<arma::mat>(C_trace[i]));
        H_trace_vec[i + 1] = (Rcpp::as<arma::mat>(H_trace[i]));
    }

    actionet::multilevel_archetypal_decomposition results =
        actionet::prune_archetypes(C_trace_vec, H_trace_vec, min_specificity_z_threshold, min_cells);

    Rcpp::List out_list;

    for (int i = 0; i < results.selected_archs.n_elem; i++)
        results.selected_archs[i]++;
    out_list["selected_archs"] = results.selected_archs;

    out_list["C_stacked"] = results.C_stacked;
    out_list["H_stacked"] = results.H_stacked;

    return out_list;
}

//' Identifies and aggregates redundant archetypes into equivalent classes
//' (Post-ACTIONet archetype processing)
//'
//' @param G Adjacency matrix of the ACTIONet graph
//' @param S_r Reduced kernel profile
//' @param archetypes Archetype profile (S*C)
//' @param C_stacked,H_stacked Output of reconstruct_archetypes()
//' @param minPoints, minClusterSize, outlier_threshold HDBSCAN parameters
//' @param reduced_dim Kernel reduction
//'
//' @return A named list: \itemize{
//' \item archetype_groups: Equivalent classes of archetypes (non-redundant)
//' \item C_unified,H_unified: C and H matrices of unified archetypes
//' \item sample_assignments: Assignment of samples/cells to unified archetypes
//' }
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
// [[Rcpp::export]]
Rcpp::List
unify_archetypes(arma::mat& S_r, arma::mat C_stacked, arma::mat H_stacked, int normalization = 0, int thread_no = 0)
{
    actionet::unification_results results =
        actionet::unify_archetypes(S_r, C_stacked, H_stacked, normalization, thread_no);

    Rcpp::List out_list;

    for (int i = 0; i < results.selected_archetypes.n_elem; i++)
        results.selected_archetypes[i]++;
    out_list["selected_archetypes"] = results.selected_archetypes;

    out_list["C_unified"] = results.C_unified;
    out_list["H_unified"] = results.H_unified;

    for (int i = 0; i < results.assigned_archetypes.n_elem; i++)
        results.assigned_archetypes[i]++;

    out_list["assigned_archetypes"] = results.assigned_archetypes;
    out_list["arch_membership_weights"] = results.arch_membership_weights;

    out_list["ontology"] = results.dag_adj;
    out_list["ontology_node_attributes"] = results.dag_node_annotations;

    return out_list;
}

//' Builds an interaction network from the multi-level archetypal decompositions
//'
//' @param H_stacked Output of the prune_archetypes() function.
//' @param density Overall density of constructed graph. The higher the density,
//' the more edges are retained (default = 1.0).
//' @param thread_no Number of parallel threads (default = 0).
//' @param mutual_edges_only Symmetrization strategy for nearest-neighbor edges.
//' If it is true, only mutual nearest-neighbors are returned (default=TRUE).
//'
//' @return G Adjacency matrix of the ACTIONet graph.
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
// [[Rcpp::export]]
arma::sp_mat buildNetwork(arma::mat H, std::string algorithm = "k*nn", std::string distance_metric = "jsd",
                          double density = 1.0, int thread_no = 0, bool mutual_edges_only = true, int k = 10)
{
    // TODO: Add as options in interface
    double M = 16, ef_construction = 200, ef = 50;
    arma::sp_mat G = actionet::buildNetwork(H, algorithm, distance_metric, density, thread_no, M, ef_construction, ef,
                                            mutual_edges_only, k);

    return G;
}

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
                         bool presmooth_network = false, double min_dist = 1, double spread = 1, double gamma = 1.0,
                         unsigned int n_epochs = 500, int thread_no = 0, int seed = 0, double learning_rate = 1.0,
                         int sim2dist = 2)
{
    arma::field<arma::mat> res =
        actionet::layoutNetwork_xmap(G, initial_position, presmooth_network, method, min_dist, spread, gamma,
                                     n_epochs, thread_no, seed, learning_rate, sim2dist);

    Rcpp::List out_list;
    out_list["coordinates"] = res(0);
    out_list["coordinates_3D"] = res(1);
    out_list["colors"] = res(2);

    return out_list;
}

//' Aggregate matrix within groups
//'
//' @param S matrix of type "dMatrix"
//' @param sample_assignments Vector of column groupings. Group labels must be continuous integers or coercible to such.
//'
//' @return S matrix with columns of values aggregated within each group of sample_assignments
//'
// [[Rcpp::export]]
arma::mat compute_grouped_rowsums(arma::sp_mat& S, arma::Col<unsigned long long> sample_assignments)
{
    arma::mat pb = actionet::compute_grouped_rowsums(S, sample_assignments);

    return pb;
}

//' Aggregate matrix within groups
//'
//' @param S matrix
//' @param sample_assignments Vector of column groupings. Group labels must be continuous integers or coercible to such.
//'
//' @return S matrix with columns of values aggregated within each group of sample_assignments
//'
// [[Rcpp::export]]
arma::mat compute_grouped_rowsums_full(arma::mat& S, arma::Col<unsigned long long> sample_assignments)
{
    arma::mat pb = actionet::compute_grouped_rowsums(S, sample_assignments);

    return pb;
}

//' Average matrix within groups
//'
//' @param S matrix of type "dMatrix"
//' @param sample_assignments Vector of column groupings. Group labels must be continuous integers or coercible to such.
//'
//' @return S matrix with columns of values average within each group of sample_assignments
//'
// [[Rcpp::export]]
arma::mat compute_grouped_rowmeans(arma::sp_mat& S, arma::Col<unsigned long long> sample_assignments)
{
    arma::mat pb = actionet::compute_grouped_rowmeans(S, sample_assignments);

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
arma::mat compute_grouped_rowmeans_full(arma::mat& S, arma::Col<unsigned long long> sample_assignments)
{
    arma::mat pb = actionet::compute_grouped_rowmeans(S, sample_assignments);

    return pb;
}

// [[Rcpp::export]]
arma::mat compute_grouped_rowvars(arma::sp_mat& S, arma::Col<unsigned long long> sample_assignments)
{
    arma::mat pb = actionet::compute_grouped_rowvars(S, sample_assignments);

    return pb;
}

// [[Rcpp::export]]
arma::mat compute_grouped_rowvars_full(arma::mat& S, arma::Col<unsigned long long> sample_assignments)
{
    arma::mat pb = actionet::compute_grouped_rowvars(S, sample_assignments);

    return pb;
}

//' Compute feature specificity (from archetype footprints)
//'
//' @param S Input matrix (sparseMatrix)
//' @param H A soft membership matrix - Typically H_unified from the unify_archetypes() function.
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//' S.norm = renormalize_input_matrix(S, cell.clusters)
//' logPvals.list = compute_archetype_feature_specificity(S.norm, unification.out$H_unified)
//' specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
Rcpp::List compute_archetype_feature_specificity(arma::sp_mat& S, arma::mat& H, int thread_no = 0)
{
    arma::field<arma::mat> res = actionet::compute_feature_specificity(S, H, thread_no);

    Rcpp::List out_list;
    out_list["archetypes"] = res(0);
    out_list["upper_significance"] = res(1);
    out_list["lower_significance"] = res(2);

    return (out_list);
}

//' Compute feature specificity (from archetype footprints)
//'
//' @param S Input matrix ("matrix" type)
//' @param H A soft membership matrix - Typically H_unified from the unify_archetypes() function.
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//' S.norm = renormalize_input_matrix(S, cell.clusters)
//' logPvals.list = compute_archetype_feature_specificity(S.norm, unification.out$H_unified)
//' specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
Rcpp::List compute_archetype_feature_specificity_full(arma::mat& S, arma::mat& H, int thread_no = 0)
{
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
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//' S.norm = renormalize_input_matrix(S, cell.clusters)
//' logPvals.list = compute_cluster_feature_specificity(S.norm, cell.clusters)
//' specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
Rcpp::List compute_cluster_feature_specificity(arma::sp_mat& S, arma::uvec sample_assignments, int thread_no = 0)
{
    arma::field<arma::mat> res = actionet::compute_feature_specificity(S, sample_assignments, thread_no);

    Rcpp::List out_list;
    out_list["average_profile"] = res(0);
    out_list["upper_significance"] = res(1);
    out_list["lower_significance"] = res(2);

    return (out_list);
}

//' Compute feature specificity (from cluster assignments)
//'
//' @param S Input matrix ("matrix")
//' @param sample_assignments Vector of cluster assignments
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//' S.norm = renormalize_input_matrix(S, cell.clusters)
//' logPvals.list = compute_cluster_feature_specificity(S.norm, cell.clusters)
//' specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
Rcpp::List compute_cluster_feature_specificity_full(arma::mat& S, arma::uvec sample_assignments, int thread_no = 0)
{
    arma::field<arma::mat> res = actionet::compute_feature_specificity(S, sample_assignments, thread_no);

    Rcpp::List out_list;
    out_list["average_profile"] = res(0);
    out_list["upper_significance"] = res(1);
    out_list["lower_significance"] = res(2);

    return (out_list);
}

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
arma::uvec compute_core_number(arma::sp_mat& G)
{
    arma::uvec core_num = actionet::compute_core_number(G);

    return (core_num);
}

//' Compute coreness of subgraph vertices induced by each archetype
//'
//' @param G Input graph
//' @param sample_assignments Archetype discretization (output of unify_archetypes())
//'
//' @return cn core-number of each graph node
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' assignments = ace$archetype.assignment
//' connectivity = compute_core_number(G, assignments)
// [[Rcpp::export]]
arma::vec compute_archetype_core_centrality(arma::sp_mat& G, arma::uvec sample_assignments)
{
    arma::vec conn = actionet::compute_archetype_core_centrality(G, sample_assignments);

    return (conn);
}

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
arma::mat compute_network_diffusion_fast(arma::sp_mat& G, arma::sp_mat& X0, int thread_no = 0, double alpha = 0.85,
                                         int max_it = 3)
{
    arma::mat Diff = actionet::compute_network_diffusion_fast(G, X0, thread_no, alpha, max_it);

    return (Diff);
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
//' specificity_scores = rowFactors(ace)[["H_unified_upper_significance"]]
//' logPvals = compute_feature_specificity(
//' specificity_scores[common.genes, ], annotations[common.genes, ]
//' )
//' rownames(logPvals) = colnames(specificity_scores)
//' colnames(logPvals) = colnames(annotations)
// [[Rcpp::export]]
Rcpp::List assess_enrichment(arma::mat& scores, arma::sp_mat& associations, int thread_no = 0)
{
    arma::field<arma::mat> res = actionet::assess_enrichment(scores, associations, thread_no);

    Rcpp::List out_list;
    out_list["logPvals"] = res(0);
    out_list["thresholds"] = res(1);

    return (out_list);
}

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
arma::mat MWM_hungarian(arma::mat& G)
{
    arma::mat G_matched = actionet::MWM_hungarian(G);

    return G_matched;
}

// [[Rcpp::export]]
arma::umat MWM_rank1(arma::vec u, arma::vec v, double u_threshold = 0, double v_threshold = 0)
{
    arma::umat pairs = actionet::MWM_rank1(u, v, u_threshold, v_threshold);

    pairs = pairs + 1;

    return (pairs);
}

// [[Rcpp::export]]
Rcpp::NumericVector run_LPA(arma::sp_mat& G, arma::vec labels, double lambda = 1, int iters = 3,
                            double sig_threshold = 3, Rcpp::Nullable<Rcpp::IntegerVector> fixed_labels_ = R_NilValue,
                            int thread_no = 0)
{
    arma::uvec fixed_labels_vec;
    if (fixed_labels_.isNotNull())
    {
        Rcpp::NumericVector fixed_labels(fixed_labels_);
        fixed_labels_vec.set_size(fixed_labels.size());
        for (int i = 0; i < fixed_labels.size(); i++)
        {
            fixed_labels_vec(i) = fixed_labels(i) - 1;
        }
    }

    arma::mat new_labels =
        actionet::LPA(G, labels, lambda, iters, sig_threshold, fixed_labels_vec, thread_no);
    // vec labels_out = arma::conv_to<arma::vec>::from(new_labels);

    return (arma2vec(new_labels));
}

// [[Rcpp::export]]
arma::mat compute_marker_aggregate_stats(arma::sp_mat& G, arma::sp_mat& S, arma::sp_mat& marker_mat,
                                         double alpha = 0.85, int max_it = 5, int thread_no = 0,
                                         bool ignore_baseline_expression = false)
{
    arma::mat stats = actionet::compute_marker_aggregate_stats(G, S, marker_mat, alpha, max_it, thread_no,
                                                               ignore_baseline_expression);

    return (stats);
}

// [[Rcpp::export]]
arma::sp_mat LSI(arma::sp_mat& X, double size_factor = 100000)
{
    arma::sp_mat TFIDF = actionet::LSI(X, size_factor);

    return (TFIDF);
}

// [[Rcpp::export]]
Rcpp::List autocorrelation_Geary(arma::sp_mat G, arma::mat scores, int normalization_method = 1, int perm_no = 30,
                                 int thread_no = 0)
{
    arma::field<arma::vec> out = actionet::autocorrelation_Geary(G, scores, normalization_method, perm_no, thread_no);

    Rcpp::List res;
    res["Geary_C"] = out[0];
    res["zscore"] = out[1];
    res["mu"] = out[2];
    res["sigma"] = out[3];

    return (res);
}

// [[Rcpp::export]]
Rcpp::List autocorrelation_Geary_full(arma::mat G, arma::mat scores, int normalization_method = 1, int perm_no = 30,
                                      int thread_no = 0)
{
    arma::field<arma::vec> out = actionet::autocorrelation_Geary(G, scores, normalization_method, perm_no, thread_no);

    Rcpp::List res;
    res["Geary_C"] = out[0];
    res["zscore"] = out[1];
    res["mu"] = out[2];
    res["sigma"] = out[3];

    return (res);
}

// [[Rcpp::export]]
Rcpp::List autocorrelation_Moran(arma::sp_mat G, arma::mat scores, int normalization_method = 1, int perm_no = 30,
                                 int thread_no = 0)
{
    arma::field<arma::vec> out = actionet::autocorrelation_Moran(G, scores, normalization_method, perm_no, thread_no);

    Rcpp::List res;
    res["Moran_I"] = out[0];
    res["zscore"] = out[1];
    res["mu"] = out[2];
    res["sigma"] = out[3];

    return (res);
}

// [[Rcpp::export]]
Rcpp::List autocorrelation_Moran_full(arma::mat G, arma::mat scores, int normalization_method = 1, int perm_no = 30,
                                      int thread_no = 0)
{
    arma::field<arma::vec> out = actionet::autocorrelation_Moran(G, scores, normalization_method, perm_no, thread_no);

    Rcpp::List res;
    res["Moran_I"] = out[0];
    res["zscore"] = out[1];
    res["mu"] = out[2];
    res["sigma"] = out[3];

    return (res);
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
arma::mat compute_network_diffusion_approx(arma::sp_mat& G, arma::mat& X0, int thread_no = 0, double alpha = 0.85,
                                           int max_it = 5, double res_threshold = 1e-8, int norm_type = 0)
{
    if (G.n_rows != X0.n_rows)
    {
        stderr_printf("Dimension mismatch: G (%dx%d) and X0 (%dx%d)\n", G.n_rows, G.n_cols, X0.n_rows, X0.n_cols);
        return (arma::mat());
    }

    arma::sp_mat P = actionet::normalize_adj(G, norm_type);
    arma::mat X = actionet::compute_network_diffusion_Chebyshev(P, X0, thread_no, alpha, max_it, res_threshold);

    return (X);
}

// [[Rcpp::export]]
arma::mat aggregate_genesets_mahalanobis_2archs(arma::sp_mat& G, arma::sp_mat& S, arma::sp_mat& marker_mat,
                                                int network_normalization_method = 0,
                                                int expression_normalization_method = 0, int gene_scaling_method = 0,
                                                double pre_alpha = 0.85, double post_alpha = 0.85, int thread_no = 0)
{
    arma::mat stats = actionet::aggregate_genesets_mahalanobis_2archs(G, S, marker_mat, network_normalization_method,
                                                                      expression_normalization_method,
                                                                      gene_scaling_method,
                                                                      pre_alpha, post_alpha, thread_no);

    return (stats);
}

// [[Rcpp::export]]
arma::mat aggregate_genesets_mahalanobis_2gmm(arma::sp_mat& G, arma::sp_mat& S, arma::sp_mat& marker_mat,
                                              int network_normalization_method = 0,
                                              int expression_normalization_method = 0, int gene_scaling_method = 0,
                                              double pre_alpha = 0.85, double post_alpha = 0.85, int thread_no = 0)
{
    arma::mat stats = actionet::aggregate_genesets_mahalanobis_2gmm(G, S, marker_mat, network_normalization_method,
                                                                    expression_normalization_method,
                                                                    gene_scaling_method,
                                                                    pre_alpha, post_alpha, thread_no);

    return (stats);
}

// TODO: Update and remove. Single reference.
// [[Rcpp::export]]
arma::mat normalize_mat(arma::mat& X, int p = 0, int dim = 0)
{
    arma::mat X_norm = actionet::normalize_matrix(X, p, dim);

    return (X_norm);
}

// TODO: Update and remove. Single reference.
// [[Rcpp::export]]
arma::sp_mat normalize_spmat(arma::sp_mat& X, int p = 0, int dim = 0)
{
    arma::sp_mat X_norm = actionet::normalize_matrix(X, p, dim);

    return (X_norm);
}

// [[Rcpp::export]]
arma::vec xicor(arma::vec xvec, arma::vec yvec, bool compute_pval = true, int seed = 0)
{
    arma::vec res = actionet::xicor(xvec, yvec, compute_pval, seed);

    return (res);
}

// [[Rcpp::export]]
Rcpp::List XICOR(arma::mat& X, arma::mat& Y, bool compute_pval = true, int seed = 0, int thread_no = 0)
{
    arma::field<arma::mat> out = actionet::XICOR(X, Y, compute_pval, seed, thread_no);

    Rcpp::List res;
    res["XI"] = out(0);
    res["Z"] = out(1);

    return (res);
}

// [[Rcpp::export]]
Rcpp::List aggregate_genesets(arma::sp_mat& G, arma::sp_mat& S, arma::sp_mat& marker_mat,
                              int network_normalization_method = 0, double alpha = 0.85, int thread_no = 0)
{
    arma::field<arma::mat> stats = actionet::aggregate_genesets_vision(G, S, marker_mat, network_normalization_method,
                                                                       alpha, thread_no);

    Rcpp::List res;

    res["stats_norm_smoothed"] = stats[0];
    res["stats_norm"] = stats[1];
    res["stats"] = stats[2];

    return (res);
}

// [[Rcpp::export]]
arma::mat assess_label_enrichment(arma::sp_mat& G, arma::mat& M, int thread_no = 0)
{
    arma::mat logPvals = actionet::assess_label_enrichment(G, M, thread_no);

    return (logPvals);
}

// [[Rcpp::export]]
Rcpp::List
autocorrelation_Moran_parametric_full(arma::mat G, arma::mat scores, int normalization_method = 4, int thread_no = 0)
{
    arma::field<arma::vec> out = actionet::autocorrelation_Moran_parametric(G, scores, normalization_method, thread_no);

    Rcpp::List res;
    res["stat"] = out[0];
    res["zscore"] = out[1];
    res["mu"] = out[2];
    res["sigma"] = out[3];

    return (res);
}

// [[Rcpp::export]]
Rcpp::List
autocorrelation_Moran_parametric(arma::sp_mat G, arma::mat scores, int normalization_method = 4, int thread_no = 0)
{
    arma::field<arma::vec> out = actionet::autocorrelation_Moran_parametric(G, scores, normalization_method, thread_no);

    Rcpp::List res;
    res["stat"] = out[0];
    res["zscore"] = out[1];
    res["mu"] = out[2];
    res["sigma"] = out[3];

    return (res);
}
