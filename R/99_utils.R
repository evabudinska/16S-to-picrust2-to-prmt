#' Small helpers and assertions

#' Ensure UTF-8 locale for reproducibility (optional)
ensure_utf8_locale <- function() {
  loc <- Sys.getlocale("LC_CTYPE")
  if (!grepl("UTF", loc, ignore.case = TRUE)) {
    message("Consider setting UTF-8 locale for consistent behavior: e.g., Sys.setlocale('LC_CTYPE','en_US.UTF-8')")
  }
}

