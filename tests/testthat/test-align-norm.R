test_that("align + relative normalization are deterministic", {
  source("R/02_io.R"); source("R/03_preprocessing.R")

  abund <- tibble::tibble(function = c("EC:1","EC:2"), S2 = c(5,5), S1 = c(5,15))
  meta  <- tibble::tibble(sample_id = c("S1","S2"), group = c("A","B"))

  aligned <- align_samples(abund, id_col = "function", metadata = meta)
  expect_equal(names(aligned$abund), c("function","S1","S2"))  # sorted
  rel <- rel_abundance(aligned$abund, id_col = "function")
  expect_equal(colSums(as.matrix(rel[, c("S1","S2")])) |> as.numeric(), c(1,1), tolerance = 1e-12)
})

