# Matrix projection to compound scores
# Refactor of your get_cmp_scores()

suppressPackageStartupMessages({
  library(data.table)
})

#' Project KO-normalized weights to compound space
#' @param emm numeric matrix/data.frame of compounds x KOs (e.g., stoichiometry or network weights)
#' @param norm_kos data.table with columns KO and sample columns (weights per KO per sample)
#' @return data.table with compounds as key and sample columns
get_cmp_scores <- function(emm, norm_kos) {
  metlen <- nrow(emm)
  nsamp  <- ncol(norm_kos) - 1L

  norm_kos_sub <- norm_kos[KO %in% colnames(emm)]
  subjects <- setdiff(colnames(norm_kos), "KO")

  # reorder EMM columns to match KOs
  emm_sub <- emm[, match(norm_kos_sub[, KO], colnames(emm)), drop = FALSE]

  if (all(colnames(emm_sub) == norm_kos_sub[, KO])) {
    cmp <- matrix(NA_real_, nrow = metlen, ncol = nsamp)
    for (m in seq_len(nsamp)) {
      cmp[, m] <- as.matrix(emm_sub) %*% as.numeric(unlist(norm_kos_sub[, subjects[m], with = FALSE]))
    }
  } else if (setequal(colnames(emm_sub), norm_kos_sub[, KO])) {
    ord <- order(colnames(emm_sub))
    emm_sub <- emm_sub[, ord, drop = FALSE]
    norm_kos_sub <- norm_kos_sub[order(KO)]
    cmp <- matrix(NA_real_, nrow = metlen, ncol = nsamp)
    for (m in seq_len(nsamp)) {
      cmp[, m] <- as.matrix(emm_sub) %*% as.numeric(unlist(norm_kos_sub[, subjects[m], with = FALSE]))
    }
  } else {
    stop("KOs in EMM and norm_kos do not match. Check inputs.")
  }

  cmp_dt <- data.table(cmp)
  setnames(cmp_dt, subjects)
  cmp_dt[, compound := rownames(emm)]
  setkey(cmp_dt, compound)
  cmp_dt
}

