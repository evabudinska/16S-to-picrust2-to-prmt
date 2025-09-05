#' Map PICRUSt2 features to PRMT sets and compute PRMT scores
#'
#' You will plug your “PICRUSt2 feature -> PRMT mapping” here.
#' Keep mapping as a long table: columns: prmt_id, feature_id, weight (optional).

#' Read PRMT mapping table
#' @param path Path to mapping (CSV/TSV) with prmt_id, feature_id, weight?
#' @return Tibble with required columns.
read_prmt_mapping <- function(path) {
  stopifnot(file.exists(path))
  map <- readr::read_delim(path, delim = "\t", col_types = readr::cols())
  req <- c("prmt_id", "feature_id")
  if (!all(req %in% names(map))) {
    stop("Mapping must contain: ", paste(req, collapse = ", "))
  }
  if (!"weight" %in% names(map)) map$weight <- 1.0
  map
}

#' Compute PRMT scores by weighted sum (default) or mean
#' @param abund_wide Tibble with id_col + sample columns (relative or raw).
#' @param id_col Feature ID column name (must match mapping$feature_id).
#' @param mapping Tibble with prmt_id, feature_id, weight.
#' @param method "sum" or "mean".
#' @return Tibble: prmt_id + sample columns (PRMT scores).
compute_prmt_scores <- function(abund_wide, id_col, mapping, method = c("sum","mean")) {
  method <- match.arg(method)
  scols  <- setdiff(names(abund_wide), id_col)

  # Join abundance with mapping to get weights per feature
  long <- abund_wide |>
    tidyr::pivot_longer(all_of(scols), names_to = "sample_id", values_to = "value") |>
    dplyr::rename(feature_id = !!id_col) |>
    dplyr::inner_join(mapping, by = "feature_id")

  if (nrow(long) == 0) {
    stop("No overlap between features in abundance and mapping.")
  }

  long <- long |>
    mutate(w_value = value * weight)

  agg <- if (method == "sum") {
    long |>
      group_by(prmt_id, sample_id) |>
      summarise(score = sum(w_value, na.rm = TRUE), .groups = "drop")
  } else {
    long |>
      group_by(prmt_id, sample_id) |>
      summarise(score = mean(w_value, na.rm = TRUE), .groups = "drop")
  }

  # Wide PRMT x samples
  prmt_wide <- tidyr::pivot_wider(agg, names_from = sample_id, values_from = score)
  prmt_wide <- arrange(prmt_wide, prmt_id)
  prmt_wide
}

