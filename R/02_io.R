#' I/O utilities for PICRUSt2 -> PRMT pipeline
#'
#' These functions read standard PICRUSt2 outputs and metadata.
#' All functions return tibbles/data.tables with well-defined columns.

#' Read a PICRUSt2 abundance table (e.g., EC or pathway predictions)
#' @param path File path to TSV/CSV.
#' @param id_col Name of the feature ID column (e.g., "function" or "EC").
#' @return A tibble with columns: id_col + sample columns (wide format).
read_picrust2_abundance <- function(path, id_col = "function") {
  stopifnot(file.exists(path))
  dt <- data.table::fread(path)
  if (!(id_col %in% names(dt))) {
    stop("`id_col` '", id_col, "' not found in file: ", path)
  }
  # Ensure unique IDs and no duplicated column names
  dt <- as_tibble(dt)
  if (anyDuplicated(dt[[id_col]]) > 0) {
    warning("Duplicated feature IDs found; keeping first occurrence.")
    dt <- dt[!duplicated(.data[[id_col]]), , drop = FALSE]
  }
  dt
}

#' Read sample metadata (sample-by-row)
#' @param path Path to metadata (CSV/TSV); must contain 'sample_id'.
#' @return Tibble with at least sample_id column.
read_metadata <- function(path) {
  stopifnot(file.exists(path))
  meta <- readr::read_delim(path, delim = guess_delim_(path), col_types = readr::cols())
  if (!"sample_id" %in% names(meta)) {
    stop("Metadata must contain a 'sample_id' column.")
  }
  meta
}

#' Guess delimiter by file extension (simple helper)
guess_delim_ <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) return(",")
  "\t"
}

