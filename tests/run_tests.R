# tests/run_tests.R — run all tests in a non-package repo
if (!requireNamespace("testthat", quietly = TRUE)) install.packages("testthat")

# (optional) set project root if spúšťaš z iného adresára
# setwd("/path/to/picrust2-to-prmt")

# Make sure testthat sees helpers
Sys.setenv("TESTTHAT" = "true")

# Source all R/ files once (namiesto source() v každom teste)
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source, chdir = TRUE))

# Run tests
testthat::test_dir("tests/testthat", reporter = "summary")

