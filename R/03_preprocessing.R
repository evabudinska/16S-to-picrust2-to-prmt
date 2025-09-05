#' Pre-processing steps:
#' - Harmonize sample IDs between abundance and metadata
#' - Optional filtering / normalization
#' - Return clean objects ready for PRMT mapping

#' Align samples between abundance and metadata
#' @param abund_tbl Tibble from read_picrust2_abundance()
#' @param id_col Feature ID column name.
#' @param metadata Tibble with sample_id column.
#' @return list(abund=filtered_wide, meta=filtered_meta)
align_samples <- function(abund_tbl, id_col, metadata) {
  sample_cols <- setdiff(names(abund_tbl), id_col)
  intersecting <- intersect(sample_cols, metadata$sample_id)
  if (!length(intersecting)) stop("No overlapping sample IDs between abundance and metadata.")
  abund_f <- abund_tbl[, c(id_col, intersecting)]
  meta_f  <- metadata %>% filter(sample_id %in% intersecting)
  # Deterministic column order
  abund_f <- abund_f[, c(id_col, sort(intersecting))]
  list(abund = abund_f, meta = arrange(meta_f, sample_id))
}

#' Minimal normalization (example): convert to relative abundances per sample
#' @param abund_wide Tibble with id_col + sample columns.
#' @param id_col Feature ID column name.
#' @return Tibble of same shape with relative abundances.
rel_abundance <- function(abund_wide, id_col) {
  scols <- setdiff(names(abund_wide), id_col)
  mat <- as.matrix(abund_wide[, scols, drop = FALSE])
  col_sums <- colSums(mat, na.rm = TRUE)
  # Avoid division by zero
  col_sums[col_sums == 0] <- 1
  rel <- sweep(mat, 2, col_sums, "/")
  out <- bind_cols(abund_wide[id_col], as_tibble(rel))
  out
}

