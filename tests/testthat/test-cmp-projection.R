test_that("get_cmp_scores projects KO weights to compounds", {
  source("R/10_prmt_cmp_scores.R")
  library(data.table)

  # EMM (compounds x KOs)
  emm <- matrix(c(
    -1, 0,   # C00001: KO1=-1, KO2=0
     1, 1    # C00002: KO1=1,  KO2=1
  ), nrow = 2, byrow = TRUE)
  rownames(emm) <- c("C00001","C00002"); colnames(emm) <- c("K1","K2")

  # norm_kos (KO + 2 samples)
  nk <- data.table::data.table(KO = c("K1","K2"), S1 = c(0.6, 0.4), S2 = c(0.0, 1.0))

  cmp <- get_cmp_scores(emm, nk)
  expect_true(all(c("S1","S2","compound") %in% names(cmp)))
  # For S1: C00001 = -1*0.6 + 0*0.4 = -0.6 ; C00002 = 1*0.6 + 1*0.4 = 1.0
  expect_equal(as.numeric(cmp[compound=="C00001","S1"]), -0.6, tolerance = 1e-12)
  expect_equal(as.numeric(cmp[compound=="C00002","S1"]),  1.0, tolerance = 1e-12)
})

