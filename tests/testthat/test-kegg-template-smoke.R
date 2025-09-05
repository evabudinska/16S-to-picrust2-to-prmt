test_that("generate_network_template_kegg parses a tiny mapformula", {
  source("R/08_network_template.R")

  # Minimal “mapformula” lines mimicking KGML output
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c(
    "R00001: 00010: C00001 <=> C00002",
    "R00002: 00020: C00002 => C00003"
  ), tmp)

  # Fake all_kegg with minimal structures to let the function run to QC point
  # We only exercise parsing/merging parts, then skip deep QC if not available.
  all_kegg <- list(
    kos_to_rxns = data.table::data.table(Rxn = c("R00001","R00002"), KO = c("K00001","K00002")),
    Reactions   = c("R00001","R00002"),
    Reaction_info = list(
      list(EQUATION = "2 C00001 <=> C00002", ORTHOLOGY = c("K00001")),
      list(EQUATION = "C00002 <=> C00003",   ORTHOLOGY = c("K00002"))
    )
  )

  tpl <- generate_network_template_kegg(tmp, all_kegg, write_out = FALSE)
  expect_true(all(c("Rxn","KO","Reac","Prod","stoichReac","stoichProd") %in% names(tpl)))
  expect_true(nrow(tpl) >= 2)
})

