test_that("compute_prmt_scores produces expected PRMT × samples matrix", {
  source("R/04_prmt_mapping.R")

  # Abundance (relative numbers for stable expectation)
  abund <- tibble::tibble(
    function = c("EC:1.1.1.1","EC:2.2.2.2","EC:9.9.9.9"),
    S1 = c(0.6, 0.4, 0.0),
    S2 = c(0.2, 0.0, 0.8)
  )
  # Mapping: PRMT_A uses EC1 (w=1) and EC2 (w=2), PRMT_B uses EC9 (w=1)
  mapping <- tibble::tibble(
    prmt_id = c("PRMT_A","PRMT_A","PRMT_B"),
    feature_id = c("EC:1.1.1.1","EC:2.2.2.2","EC:9.9.9.9"),
    weight = c(1,2,1)
  )

  prmt <- compute_prmt_scores(abund, id_col = "function", mapping = mapping, method = "sum")
  # Expected: PRMT_A,S1 = 1*0.6 + 2*0.4 = 1.4 ; PRMT_A,S2 = 1*0.2 + 2*0.0 = 0.2
  #           PRMT_B,S1 = 0 ; PRMT_B,S2 = 0.8
  expect_equal(prmt$S1[prmt$prmt_id=="PRMT_A"], 1.4, tolerance = 1e-12)
  expect_equal(prmt$S2[prmt$prmt_id=="PRMT_A"], 0.2, tolerance = 1e-12)
  expect_equal(prmt$S1[prmt$prmt_id=="PRMT_B"], 0.0, tolerance = 1e-12)
  expect_equal(prmt$S2[prmt$prmt_id=="PRMT_B"], 0.8, tolerance = 1e-12)
})

