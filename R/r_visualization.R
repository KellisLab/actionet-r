#' @export
layoutNetwork <- function(
    adata = NULL,
    initial_coordinates = NULL,
    net_slot = "actionet",
    layer = "logcounts",
    method = c("umap", "tumap", "largevis", "leopold", "leopold2"),
    n_components = 2,
    spread = 1.0,
    min_dist = 1.0,
    n_epochs = NULL,
    learning_rate = 1.0,
    repulsion_strength = 1.0,
    negative_sample_rate = 3.0,
    approx_pow = TRUE,
    pcg_rand = TRUE,
    rng_type = NULL,
    batch = TRUE,
    grain_size = 1,
    ai = NULL,
    aj = NULL,
    seed = 0,
    thread_no = 0,
    verbose = TRUE,
    a = NULL,
    b = NULL,
    opt_method = c("adam", "sgd"),
    alpha = learning_rate,
    beta1 = 0.5,
    beta2 = 0.9,
    eps = 1e-7,
    repair_disconnected = TRUE,
    map_slot_out = NULL,
    return_raw = FALSE,
    assay_name = NULL,
    obj = NULL) {
    adata <- .resolve_container_arg(adata = adata, obj = obj)
    layer <- .resolve_layer_arg(
        layer = layer,
        assay_name = assay_name,
        default = "logcounts",
        layer_missing = missing(layer),
        assay_name_missing = missing(assay_name)
    )
    force(alpha)
    method <- tolower(method)
    method <- match.arg(method, several.ok = TRUE)[1]

    opt_method <- tolower(opt_method)
    opt_method <- match.arg(opt_method, several.ok = TRUE)[1]

    if (!is.null(rng_type)) {
        if (!is.character(rng_type) || length(rng_type) != 1) {
            stop("'rng_type' must be a single character value")
        }
        rng_type <- tolower(rng_type)
        if (!rng_type %in% c("pcg", "tausworthe", "deterministic")) {
            stop("'rng_type' must be one of: 'pcg', 'tausworthe', 'deterministic'")
        }
    }

    is_ace <- .validate_ace(adata, allow_se_like = TRUE, error_on_fail = FALSE, return_elem = FALSE)
    G <- .ace_or_net(
        obj = adata,
        net_slot = net_slot,
        matrix_type = "sparse",
        sparse_type = "CsparseMatrix",
        force_type = TRUE,
        obj_name = "ace"
    )

    if (!is.null(initial_coordinates)) {
        if (is.matrix(initial_coordinates)) {
            initial_coordinates <- .validate_matrix(
                x = initial_coordinates,
                matrix_type = "dense",
                force_type = TRUE,
                return_elem = TRUE
            )
        } else if (is_ace) {
            initial_coordinates <- .validate_map(
                adata,
                map_slot = initial_coordinates,
                matrix_type = "dense",
                force_type = TRUE,
                return_elem = TRUE
            )
        } else {
            err <- sprintf("'initial_coordinates' must be type 'matrix' for 'obj' type '%s'", class(adata))
            stop(err)
        }
    } else {
        if (!is_ace) {
            err <- sprintf("'initial_coordinates' cannot be NULL for 'obj' type '%s'", class(adata))
            stop(err)
        } else {
            msg <- sprintf("Computing initial coordinates from assay '%s'", layer)
            message(msg)
            svd.out <- runSVD(
                X = .validate_assay(adata, assay_name = layer, return_elem = TRUE),
                k = base::max(3, n_components),
                seed = seed,
                verbose = verbose
            )
            # Post-flip: assay is cells x genes, so u is the cell-space (obs x k) embedding.
            initial_coordinates <- scale(svd.out$u)
        }
    }

    if (NROW(initial_coordinates) != .n_obs(adata)) {
        err <- sprintf("'NROW(initial_coordinates)' (%d) does not match number of cells in object (%d)", NROW(initial_coordinates), .n_obs(adata))
        stop(err)
    }

    if (n_components < 2) {
        err <- sprintf("'n_components' (%d) must be >= 2", n_components)
        stop(err)
    }

    if (NCOL(initial_coordinates) < n_components) {
        err <- sprintf("'NCOL(initial_coordinates)' must be >= 'n_components' (%d)", n_components)
        stop(err)
    }

    .validate_optional_coeff <- function(values, name) {
        if (is.null(values)) {
            return(NULL)
        }
        if (!is.numeric(values)) {
            err <- sprintf("'%s' must be numeric when provided", name)
            stop(err)
        }
        values <- as.numeric(values)
        if (length(values) != .n_obs(adata)) {
            err <- sprintf("'%s' must have length %d (number of observations)", name, .n_obs(adata))
            stop(err)
        }
        values
    }

    ai <- .validate_optional_coeff(ai, "ai")
    aj <- .validate_optional_coeff(aj, "aj")

    if (method == "leopold" && is.null(ai)) {
        stop("'ai' must be provided when method='leopold'")
    }
    if (method == "leopold2" && (is.null(ai) || is.null(aj))) {
        stop("'ai' and 'aj' must be provided when method='leopold2'")
    }

    if (is.null(a) || is.null(b)) {
        a <- b <- 0
    }

    embedding <- C_layoutNetwork(
        G = G,
        initial_coordinates = initial_coordinates,
        method = method,
        n_components = n_components,
        spread = spread,
        min_dist = min_dist,
        n_epochs = ifelse(is.null(n_epochs), 0, n_epochs),
        learning_rate = learning_rate,
        repulsion_strength = repulsion_strength,
        negative_sample_rate = negative_sample_rate,
        approx_pow = approx_pow,
        pcg_rand = pcg_rand,
        rng_type = if (is.null(rng_type)) "" else rng_type,
        batch = batch,
        grain_size = grain_size,
        ai = if (is.null(ai)) numeric(0) else ai,
        aj = if (is.null(aj)) numeric(0) else aj,
        seed = seed,
        thread_no = thread_no,
        verbose = verbose,
        a = a,
        b = b,
        opt_method = opt_method,
        alpha = alpha,
        beta1 = beta1,
        beta2 = beta2,
        eps = eps,
        repair_disconnected = repair_disconnected
    )

    if (is_ace && !return_raw) {
        rownames(embedding) <- .actionet_colnames(adata)
        if (is.null(map_slot_out)) {
            map_slot_out <- sprintf("%s_%dd_%s", method, n_components, net_slot)
        }
        colMaps(adata)[[map_slot_out]] <- embedding
        colMapTypes(adata)[[map_slot_out]] <- "embedding"
        return(adata)
    }
    return(embedding)
}

#' @export
computeNodeColors <- function(
    adata = NULL,
    embedding_slot = "umap_3d_actionet",
    color_slot_out = NULL,
    thread_no = 1,
    return_raw = FALSE,
    obj = NULL) {
    adata <- .resolve_container_arg(adata = adata, obj = obj)
    is_ace <- .validate_ace(adata, allow_se_like = TRUE, error_on_fail = FALSE, return_elem = FALSE)

    coordinates <- .ace_or_map(
        obj = adata,
        map_slot = embedding_slot,
        matrix_type = "dense",
        force_type = TRUE,
        transpose_map = FALSE,
        return_elem = TRUE
    )

    colors <- C_computeNodeColors(coordinates, thread_no)

    if (is_ace && !return_raw) {
        if (is.null(color_slot_out)) {
            color_slot_out <- sprintf("colors_%s", embedding_slot)
        }
        colMaps(adata)[[color_slot_out]] <- colors
        return(adata)
    }
    return(colors)
}
