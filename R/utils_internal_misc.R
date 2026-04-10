.get_features <- function(obj, features_use = NULL, allow_empty = FALSE, features_name = "features_use") {
    if (is.null(features_use) || !.is_se_like(obj)) {
        features_use <- if (.is_anndata(obj)) .actionet_rownames(obj) else rownames(obj)
    } else {
        features_use <- .validate_vector_attr(
            obj,
            attr = features_use,
            return_type = "data",
            dim = 1,
            return_elem = TRUE
        )
    }

    if (length(features_use) == 0) {
        if (!allow_empty) {
            err <- sprintf("%s cannot be empty", features_name)
            stop(err)
        }
    }
    return(features_use)
}

## Deprecated
.get_feature_vec <- function(ace, features_use = NULL) {
    if (is.null(features_use)) {
        features_use <- .actionet_rownames(ace)
    } else {
        features_use <- .validate_vector_attr(
            obj = ace,
            attr = features_use,
            return_type = "data",
            dim = 1
        )
    }
    return(features_use)
}

.preprocess_annotation_labels <- function(labels, ace = NULL) {
    if (is.null(labels)) {
        return(NULL)
    }

    if (is.character(labels)) {
        labels <- factor(.validate_vector_attr(ace, attr = labels, return_type = "data"))
    }

    if ((length(labels) > 1) & is.logical(labels)) {
        labels <- factor(as.numeric(labels), levels = c(0, 1), labels = c("No", "Yes"))
    }

    if (is.factor(labels)) {
        v <- as.numeric(labels)
        names(v) <- levels(labels)[v]
        labels <- v
    }
    if (is.matrix(labels)) {
        L <- as.numeric(labels)
        names(L) <- names(labels)
        labels <- L
    }

    if (is.null(names(labels)) | length(unique(names(labels))) > length(labels) / 2) {
        names(labels) <- as.character(labels)
    }

    return(labels)
}

.preprocess_design_matrix_and_var_names <- function(design_mat, variable_name = NULL) {
    if (is.null(variable_name)) {
        variable_name <- colnames(design_mat)[ncol(design_mat)]
    }
    vn_idx <- which(variable_name == colnames(design_mat))[1]
    colnames(design_mat) <- make.names(colnames(design_mat), unique = TRUE, allow_ = FALSE)
    variable_name <- colnames(design_mat)[vn_idx]

    out <- list(design_mat = design_mat, variable_name = variable_name)
    return(out)
}
