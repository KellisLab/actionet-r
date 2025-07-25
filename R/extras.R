# HGT tail bound
Kappa <- function(p, q) {

    kl = array(1, length(p))

    suppressWarnings({ a = p * log(p/q) })
    a[p == 0] = 0

    suppressWarnings({ b = (1 - p) * log((1 - p)/(1 - q)) })
    b[p == 1] = 0
    k = a + b

    return(k)
}


HGT_tail <- function(
  population.size,
  success.count,
  sample.size,
  observed.success
) {

    if (sum(success.count) == 0)
        return(rep(0, length(success.count)))

    success.rate = success.count/population.size
    expected.success = sample.size * success.rate
    delta = (observed.success/expected.success) - 1

    log.tail_bound = sample.size * Kappa((1 + delta) * success.rate, success.rate)
    log.tail_bound[delta < 0] = 0
    log.tail_bound[is.na(log.tail_bound)] = 0

    return(log.tail_bound)
}


combine.logPvals <- function(
  logPvals,
  top.len = NULL,
  base = 10
) {

    if (is.null(top.len)) {
        top.len = nrow(logPvals)
    }
    kappa_val = 1/log(exp(1), base = base)
    logPvals = kappa_val * logPvals

    combined.log.pvals = -apply(logPvals, 2, function(lx) {
        perm = order(lx, decreasing = TRUE)

        return(log(top.len) - matrixStats::logSumExp(lx[perm[1:top.len]]))
    })

    return(combined.log.pvals)
}


doubleNorm <- function(
  Enrichment,
  log.transform = TRUE,
  min.threshold = 0
) {

    Enrichment[Enrichment < min.threshold] = 0

    if ((max(Enrichment) > 100) & (log.transform == TRUE)) {
        Enrichment = log1p(Enrichment)
    }
    Enrichment[is.na(Enrichment)] = 0

    rs = sqrt(Matrix::rowSums(Enrichment))
    rs[rs == 0] = 1
    D_r = Matrix::Diagonal(nrow(Enrichment), 1/rs)

    cs = sqrt(Matrix::colSums(Enrichment))
    cs[cs == 0] = 1
    D_c = Matrix::Diagonal(ncol(Enrichment), 1/cs)

    Enrichment.scaled = as.matrix(D_r %*% Enrichment %*% D_c)

    Enrichment.scaled = Enrichment.scaled/max(Enrichment.scaled)

    return(Enrichment.scaled)
}
