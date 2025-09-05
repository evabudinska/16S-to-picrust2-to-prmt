#' Optional downstream stats on PRMT scores

#' Z-score PRMT scores across samples (per PRMT)
#' @param prmt_wide Tibble with prmt_id + sample columns
#' @return Tibble with same shape (z-scored per row).
zscore_prmt <- function(prmt_wide) {
  stopifnot("prmt_id" %in% names(prmt_wide))
  scols <- setdiff(names(prmt_wide), "prmt_id")
  mat <- as.matrix(prmt_wide[, scols, drop = FALSE])
  mu  <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, sd, na.rm = TRUE)
  sdv[sdv == 0] <- 1
  z   <- (mat - mu) / sdv
  bind_cols(prmt_id = prmt_wide$prmt_id, as_tibble(z))
}

