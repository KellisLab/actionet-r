.get_mito_feats <- function(
    id_type = c("gene_name", "ensembl_id"),
    species = c("hsapiens", "mmusculus", "human", "mouse"),
    protein_coding = FALSE) {
  id_type <- tolower(id_type)
  id_type <- match.arg(id_type, several.ok = TRUE)[1]

  species <- tolower(species)
  species <- match.arg(species, several.ok = TRUE)[1]

  if (species == "human") {
    species <- "hsapiens"
  }

  if (species == "mouse") {
    species <- "mmusculus"
  }

  feats_tbl <- int_feats_mito[int_feats_mito[["species"]] == species, , drop = FALSE]

  if (protein_coding) {
    feats_tbl <- feats_tbl[feats_tbl[["feature_type"]] == "protein_coding", , drop = FALSE]
  }

  out <- feats_tbl[[id_type]]
  return(out)
}

#' @export
getFeatureAbundance <- function(
    obj,
    features,
    by = NULL,
    groups_use = NULL,
    features_use = NULL,
    metric = c("fraction", "percent", "ratio", "counts"),
    assay_name = "counts") {
  metric <- match.arg(metric, several.ok = TRUE)[1]

  is_ace <- .validate_ace(obj, allow_se_like = TRUE, error_on_fail = FALSE, return_elem = FALSE, obj_name = "obj")
  X <- .ace_or_assay(obj, assay_name = assay_name, allow_se_like = TRUE, return_elem = TRUE)

  if (is_ace) {
    features_use <- .get_feature_vec(obj, features_use = features_use)
  } else {
    if (is.null(features_use)) {
      features_use <- rownames(obj)
    } else {
      if (length(features_use) != NROW(obj)) {
        err <- sprintf("length(features_use) (%d) does not match NROW(obj) (%d)", length(features_use), NROW(obj))
        stop(err)
      }
    }
  }

  idx_feat <- which(features_use %in% features)
  if (length(idx_feat) < 1) {
    stop("No matching features")
  }

  cs_mat <- Matrix::colSums(X)
  cs_mm <- C_fastSpMatViewSum(
    X = X,
    idx = idx_feat,
    dim = 2,
    threads = .get_num_threads(base::ceiling(length(idx_feat) / 10), 0)
  )

  if (!is.null(by)) {
    if (is_ace) {
      by.split <- ACTIONetExperiment::get.data.or.split(
        obj,
        attr = by,
        groups_use = groups_use,
        to_return = "split"
      )
    } else {
      if (length(by) != NCOL(obj)) {
        err <- sprintf("length(by) (%d) does not match NCOL(obj) (%d)", length(by), NCOL(obj))
        stop(err)
      }
      by.split <- split(seq_len(NCOL(X)), f = by, drop = TRUE)
    }

    if (metric %in% c("fraction", "percent")) {
      out <- lapply(by.split, function(idx) {
        x <- cs_mm[idx] / cs_mat[idx]
        if (metric == "percent") {
          x <- 100 * x
        }
        return(x)
      })
    } else if (metric == "ratio") {
      out <- lapply(by.split, function(idx) {
        x <- cs_mm[idx] / (cs_mat[idx] - cs_mm[idx])
      })
    } else {
      out <- lapply(by.split, function(idx) {
        x <- cs_mm[idx]
      })
    }
  } else {
    if (metric %in% c("fraction", "percent")) {
      out <- cs_mm / cs_mat
      if (metric == "percent") {
        out <- 100 * out
      }
    } else if (metric == "ratio") {
      out <- cs_mm / (cs_mat - cs_mm)
    } else {
      out <- cs_mm
    }
    return(out)
  }
}


#' @export
getMitoAbundance <- function(
    obj,
    by = NULL,
    groups_use = NULL,
    features_use = NULL,
    id_type = c("gene_name", "ensembl_id"),
    species = c("hsapiens", "mmusculus", "human", "mouse"),
    protein_coding = FALSE,
    metric = c("fraction", "percent", "ratio", "counts"),
    assay_name = "counts") {
  feats_mito <- .get_mito_feats(id_type = id_type, species = species, protein_coding = protein_coding)

  out <- getFeatureAbundance(
    obj,
    features = feats_mito,
    by = by,
    groups_use = groups_use,
    features_use = features_use,
    metric = metric
  )
  return(out)
}

#' @export
plot.mtRNA.dist.by.attr <- function(
    ace,
    by,
    groups_use = NULL,
    features_use = NULL,
    assay = "counts",
    log_scale = FALSE,
    species = c("hsapiens", "mmusculus"),
    metric = c("pct", "ratio", "counts"),
    to_return = c("plot", "data"),
    palette = NULL,
    x_label = NULL,
    y_label = NULL,
    plot_title = NULL) {
  require(stats)
  to_return <- match.arg(to_return)

  frac.list <- computeMitoStats(
    ace = ace,
    by = by,
    groups_use = groups_use,
    features_use = features_use,
    assay = assay,
    species = species,
    metric = metric
  )

  df <- lapply(1:length(frac.list), function(l) data.frame(attr = names(frac.list)[l], frac = frac.list[[l]]))
  df <- do.call(rbind, df)

  if (metric == "pct") {
    y_label <- ifelse(is.null(y_label), ylab("% mtRNA"), y_label)
  } else if (metric == "ratio") {
    y_label <- ifelse(is.null(y_label), ylab("mtRNA:non-mtRNA"), y_label)
  } else {
    if (log_scale) {
      df$frac[df$frac == 0] <- 1
      df$frac <- log10(df$frac)
      y_label <- ifelse(is.null(y_label), ylab("Fraction mtRNA (log10)"), y_label)
    } else {
      y_label <- ifelse(is.null(y_label), ylab("Fraction mtRNA"), y_label)
    }
  }


  if (to_return == "data") {
    return(df)
  } else {
    p <- .plot_gg_violin(
      df,
      x = df$attr,
      y = df$frac,
      fill = df$attr,
      x_label = x_label,
      y_label = y_label,
      plot_title = plot_title,
      palette = palette
    )
    return(p)
  }
}

#' @export
plot.counts.by.attr <- function(
    ace,
    attr,
    nonzero = FALSE,
    log_scale = FALSE,
    assay = "counts",
    to_return = c("plot", "data"),
    palette = NULL,
    x_label = NULL,
    y_label = NULL,
    plot_title = NULL) {
  to_return <- match.arg(to_return)

  require(ggplot2)

  IDX <- ACTIONetExperiment::get.data.or.split(ace, attr = attr, to_return = "split")
  mat <- assays(ace)[[assay]]

  if (nonzero == TRUE) {
    mat <- mat > 0
  }

  cs_mat <- Matrix::colSums(mat)
  sums.list <- lapply(IDX, function(idx) {
    cs_mat[idx]
  })

  df <- lapply(1:length(sums.list), function(l) data.frame(attr = names(sums.list)[l], sums = sums.list[[l]]))
  df <- do.call(rbind, df)

  if (log_scale) {
    df$sums[df$sums == 0] <- 1
    df$sums <- log10(df$sums)
    y_label <- ifelse(is.null(y_label), ylab("sums (log10)"), y_label)
  } else {
    y_label <- ifelse(is.null(y_label), ylab("sums"), y_label)
  }

  if (to_return == "data") {
    return(df)
  } else {
    p <- .plot_gg_violin(
      df,
      x = df$attr,
      y = df$sums,
      fill = df$attr,
      x_label = x_label,
      y_label = y_label,
      plot_title = plot_title,
      palette = palette
    )
    return(p)
  }
}

#' @import ggplot2
.plot_gg_violin <- function(df, x, y, fill, x_label = NULL, y_label = NULL, plot_title = NULL, palette = NULL) {
  require(ggplot2)
  if (is.null(x_label)) {
    x_label <- element_blank()
  }
  if (is.null(y_label)) {
    y_label <- element_blank()
  }
  if (is.null(palette)) {
    palette <- CPal_default[seq.int(unique(x))]
  }

  all_labels <- unique(x) %>% sort()

  p <- ggplot(df, aes(x = x, y = y, fill = fill)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = "width") +
    scale_fill_manual(values = palette) +
    scale_x_discrete(
      labels = all_labels,
      position = "bottom"
    ) +
    labs(x = x_label, y = y_label) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "none"
    )

  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
  }

  return(p)
}
