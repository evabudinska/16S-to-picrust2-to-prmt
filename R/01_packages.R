#' Load and check required packages
#'
#' Keep dependencies minimal and explicit.
#' Called by pipeline/run_pipeline.R before anything else.
suppressPackageStartupMessages({
  pkgs <- c(
    "data.table",  # fast I/O and frames
    "readr",       # robust CSV/TSV reading
    "dplyr",       # data manipulation
    "tidyr",       # reshape
    "stringr",     # string utilities
    "ggplot2",     # plots
    "yaml",        # read config
    "purrr"        # functionals
    # "ComplexHeatmap", # optional; comment-in when used
    # "limma"          # optional; batch correction
  )
})

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       ". Please install them before running the pipeline.")
}

# Attach only what you need (avoid masking):
invisible(lapply(c("data.table", "dplyr", "tidyr", "stringr", "ggplot2", "yaml", "purrr"), library, character.only = TRUE))

options(stringsAsFactors = FALSE)

