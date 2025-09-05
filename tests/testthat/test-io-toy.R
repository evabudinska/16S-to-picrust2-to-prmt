test_that("PICRUSt2 readers handle EC/KO vs pathways id_col", {
  source("R/02_io.R")

  # Simulate KO/EC table shape (first column is 'function')
  ec_tbl <- tibble::tibble(function = c("EC:1.1.1.1", "EC:2.2.2.2"), S1 = c(10, 0), S2 = c(0, 5))
  tf <- tempfile(fileext = ".tsv"); readr::write_tsv(ec_tbl, tf)
  got <- read_picrust2_abundance(tf, id_col = "function")
  expect_has_cols(got, c("function", "S1", "S2"))

  # Simulate pathways shape (first column is 'pathway')
  path_tbl <- tibble::tibble(pathway = c("PWY-001", "PWY-002"), S1 = c(1, 2), S2 = c(3, 4))
  tf2 <- tempfile(fileext = ".tsv"); readr::write_tsv(path_tbl, tf2)
  got2 <- read_picrust2_abundance(tf2, id_col = "pathway")
  expect_has_cols(got2, c("pathway", "S1", "S2"))
})

