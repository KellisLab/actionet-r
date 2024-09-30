.create_design_formula <- function(vars) {
  fml <- Reduce(function(x, y) paste(x, y, sep = " + "), vars)
  fml <- paste("~", fml)
  fml <- as.formula(fml)
  return(fml)
}

.make_design_mat <- function(design, data = NULL) {
  if (is(design, "formula")) {
    if (is.null(data)) {
      stop("'data' is missing.")
    }
    design_mat <- stats::model.matrix(
      object = terms(design, keep.order = TRUE),
      data = data
    )
  } else if (is.matrix(design)) {
    design_mat <- design
  }
  return(design_mat)
}
