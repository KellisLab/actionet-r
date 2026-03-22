#' Gene expression imputation using network diffusion.
#'
#' @param adata AnnData or compatible container containing output of `runACTIONet()`.
#' @param features The list of features to perform imputation on.
#' @param features_use A vector of features of length `n_vars(adata)` or the name of a column of `adata$var` containing the features given in `features`.
#' @param alpha Depth of diffusion between (0, 1).
#' The larger it is, the deeper the diffusion, which results in less nonzeros (default = 0.85).
#' @param thread_no Number of parallel threads
#' @param diffusion_it Number of diffusion iterations (default = 5)
#' @param layer Layer in the container with normalized counts. `NULL` uses `adata$X`.
#'
#' @return Imputed feature expression matrix. Column names are set with imputed feature names and rows are observations.
#'
#' @examples
#' adata <- runACTIONet(adata)
#' imputed.features <- imputeFeatures(adata, c("CD14", "CD19", "CD3G"))

#' @export
imputeFeatures <- function(
    adata = NULL,
    features,
    algorithm = c("actionet", "pca"),
    features_use = NULL,
    alpha = 0.85,
    norm_method = "pagerank_sym",
    thread_no = 0,
    max_it = 5,
    layer = "logcounts",
    reduction_slot = "action",
    net_slot = "actionet",
    assay_name = NULL,
    ace = NULL) {
    algorithm <- match.arg(algorithm)
    adata <- .resolve_container_arg(adata = adata, ace = ace)
    layer <- .resolve_layer_arg(
      layer = layer,
      assay_name = assay_name,
      default = "logcounts",
      layer_missing = missing(layer),
      assay_name_missing = missing(assay_name)
    )

    features_use <- .get_features(adata, features_use = features_use, allow_empty = FALSE)
    matched_feat <- intersect(unique(features), features_use)
    idx_feat <- match(matched_feat, features_use)

    if (length(idx_feat) == 0) {
        err <- sprintf("No 'features' in 'features_use'")
        stop(err)
    }

    adata <- .validate_ace(adata, allow_se_like = FALSE, allow_null = FALSE, obj_name = "adata", as_ace = TRUE, return_elem = TRUE, error_on_fail = TRUE)

    X0 <- .ace_or_assay(
        adata,
        assay_name = layer,
        allow_se_like = FALSE,
        return_elem = TRUE
    )[, idx_feat, drop = FALSE]


    if (algorithm == "pca") {
        pc_smooth <- smoothKernel(
            adata = adata,
            norm_method = norm_method,
            alpha = alpha,
            max_it = max_it,
            reduction_slot = reduction_slot,
            net_slot = net_slot,
            thread_no = thread_no,
            return_raw = TRUE
        )

        H <- pc_smooth$H
        W <- pc_smooth$SVD.out$u
        W <- W[idx_feat, , drop = FALSE]

        out <- Matrix::t(W %*% Matrix::t(H))  # cells x features
        out[out < 0] <- 0
        # } else if (algorithm == "action") { # TODO: Fix this!! We need to also impute C. What alpha values?
        #     if (!("archetype_footprint" %in% names(colMaps(ace))) | (force_reimpute == TRUE)) {
        #         H <- networkDiffusion(
        #             obj = ace,
        #             scores = colMaps(ace)[["H_merged"]],
        #             norm_method = norm_method,
        #             alpha = alpha,
        #             thread_no = thread_no,
        #             max_it = max_it,
        #             net_slot = net_slot
        #         )
        #     } else {
        #         H <- ace$archetype_footprint
        #     }
        #     C <- colMaps(ace)$C_merged
        #     W <- as.matrix(X0 %*% C)
        #     out <- W %*% Matrix::t(H)
    } else {
        out <- networkDiffusion(
            adata = adata,
            scores = X0,
            norm_method = norm_method,
            alpha = alpha,
            thread_no = thread_no,
            approx = TRUE,
            max_it = max_it,
            net_slot = net_slot,
            return_raw = TRUE
        )
    }

    # Re-scale expression of features (out and X0 are both cells x features)
    m1 <- apply(X0, 2, max)
    m2 <- apply(out, 2, max)
    ratio <- m1 / m2
    ratio[m2 == 0] <- 1
    D <- Matrix::Diagonal(NCOL(out), ratio)
    out <- as.matrix(out %*% D)

    colnames(out) <- matched_feat
    rownames(out) <- .actionet_colnames(adata)
    return(out)
}


#' Imputing expression of genes by interpolating over archetype profile
#'
#' @param adata AnnData or compatible ACTIONet output container.
#' @param genes List of genes to impute
#' @param features_use A vector of features of length `n_vars(adata)` or the name of a column of `adata$var` containing the genes given in `genes`.
#'
#' @return A matrix of imputed expression values
#'
#' @examples
#' expression_imputed <- impute.genes.using.archetypes(adata, genes)
#' @export
impute.genes.using.archetypes <- function(adata = NULL, genes, features_use = NULL, ace = NULL) {
    adata <- .resolve_container_arg(adata = adata, ace = ace)
    adata <- .validate_ace(adata, as_ace = TRUE, allow_se_like = TRUE, return_elem = TRUE)
    features_use <- .get_feature_vec(adata, features_use = features_use)
    matched_feat <- intersect(unique(genes), features_use)
    idx_feat <- match(matched_feat, features_use)

    Z <- rowMaps(adata)[["archetype_feat_profile"]][idx_feat, , drop = FALSE]
    H <- Matrix::t(colMaps(adata)[["H_merged"]])  # cells x archetypes → archetypes x cells for Z %*% H

    expression_imputed <- Matrix::t(Z %*% H)
    colnames(expression_imputed) <- matched_feat
    rownames(expression_imputed) <- .actionet_colnames(adata)

    return(expression_imputed)
}


#' Imputing expression specificity of genes by interpolating over archetype profile
#'
#' @param adata AnnData or compatible ACTIONet output container.
#' @param genes List of genes to impute
#' @param features_use A vector of features of length `n_vars(adata)` or the name of a column of `adata$var` containing the genes given in `genes`.
#'
#' @return A matrix of imputed expression values
#'
#' @examples
#' expression_imputed <- impute.specific.genes.using.archetypes(adata, genes)
#' @export
impute.specific.genes.using.archetypes <- function(adata = NULL, genes, features_use = NULL, ace = NULL) {
    adata <- .resolve_container_arg(adata = adata, ace = ace)
    adata <- .validate_ace(adata, as_ace = TRUE, allow_se_like = TRUE, return_elem = TRUE)
    features_use <- .get_feature_vec(adata, features_use = features_use)
    matched_feat <- intersect(unique(genes), features_use)
    idx_feat <- match(matched_feat, features_use)

    Z <- log1p(rowMaps(adata)[["archetype_feat_specificity_upper"]][idx_feat, , drop = FALSE])
    H <- Matrix::t(colMaps(adata)[["H_merged"]])  # cells x archetypes → archetypes x cells for Z %*% H

    expression_imputed <- Matrix::t(Z %*% H)
    colnames(expression_imputed) <- matched_feat
    rownames(expression_imputed) <- .actionet_colnames(adata)

    return(expression_imputed)
}
