# Build KEGG-based reaction template expanded to (Rxn, KO, Reac, Prod, stoichiometries)
# Derived from your generate_network_template_kegg(), refactored and commented.

suppressPackageStartupMessages({
  library(data.table)
})

#' Generate network template from reaction_mapformula & all_kegg
#' @param mapformula_file colon-separated "Rxn: Path: Reac <=> Prod"
#' @param all_kegg list returned by get_kegg_reaction_info()
#' @param write_out logical; write TSV with timestamped filename
#' @return data.table with columns: Rxn, Path, ReacProd, Reac, Prod, KO, stoichReac, stoichProd
generate_network_template_kegg <- function(mapformula_file, all_kegg, write_out = TRUE) {
  # Read and pre-process reaction mapformula
  mf <- fread(mapformula_file, colClasses = "character", sep = ":",
              header = FALSE, quote = "", fill = TRUE)
  setnames(mf, c("Rxn", "Path", "ReacProd"))
  mf[, Rxn := trimws(Rxn)]
  mf[, Path := trimws(Path)]
  mf[, ReacProd := trimws(ReacProd)]

  # Split reactants / products; keep letters-only tokens (compound IDs)
  mf[, Reac := lapply(ReacProd, function(x) unlist(strsplit(unlist(strsplit(x, "="))[1], " ")))]
  mf[, Prod := lapply(ReacProd, function(x) unlist(strsplit(unlist(strsplit(x, "="))[2], " ")))]
  mf[, Reac := sapply(Reac, function(x) x[grepl("[A-Za-z]", x)])]
  mf[, Prod := sapply(Prod, function(x) x[grepl("[A-Za-z]", x)])]

  # Duplicate reversed direction for reversible reactions; normalize backward arrows
  for (j in seq_len(nrow(mf))) {
    submap <- mf[j]
    if (grepl("<=>", submap[, ReacProd])) {
      setnames(submap, c("Reac", "Prod"), c("Prod", "Reac"))
      submap <- submap[, .(Rxn, Path, ReacProd, Reac, Prod)]
      mf <- rbind(mf, submap)
    }
    if (grepl("<= ", submap[, ReacProd])) {
      setnames(submap, c("Reac", "Prod"), c("Prod", "Reac"))
      submap <- submap[, .(Rxn, Path, ReacProd, Reac, Prod)]
      mf[j] <- submap
    }
  }

  # Expand multiple reactants/products into all combinations
  new_mf <- data.table(Rxn = character(), Path = character(),
                       ReacProd = character(), Reac = character(), Prod = character())
  for (j in seq_len(nrow(mf))) {
    reacs <- mf[j, Reac][[1]]
    prods <- mf[j, Prod][[1]]
    submap <- mf[j]
    for (k in seq_along(reacs)) {
      for (m in seq_along(prods)) {
        submap[, Reac := reacs[k]]
        submap[, Prod := prods[m]]
        new_mf <- rbind(new_mf, submap)
      }
    }
  }
  mf <- new_mf

  # Filter general map 01100
  mf <- mf[Path != " 01100" & Path != 1100 & Path != "01100"]
  # Link reactions to KOs
  mf <- merge(mf, all_kegg$kos_to_rxns, by = "Rxn", all.x = TRUE, allow.cartesian = TRUE)
  mf <- mf[!is.na(KO)]
  setkey(mf, NULL)

  rxn_table <- unique(mf)
  rxn_check <- unique(rxn_table[, .(Rxn, KO, Reac, Prod)])

  # QC against KEGG reaction annotation: KO present & both compounds included in EQUATION
  rxn_id_check <- sapply(seq_len(nrow(rxn_check)), function(i) {
    ind <- which(all_kegg$Reactions == rxn_check[i, Rxn])
    if (length(ind) < 1) {
      return("noMatch")
    } else {
      info <- all_kegg$Reaction_info[[ind]]
      has_ko <- (rxn_check[i, KO] %in% names(info$ORTHOLOGY)) ||
                (rxn_check[i, KO] %in% info$ORTHOLOGY) ||
                (rxn_check[i, KO] %in% gsub("KO: ", "", info$ORTHOLOGY))
      has_cmpds <- (!is.null(info$EQUATION)) &&
        any(grepl(rxn_check[i, Reac], info$EQUATION)) &&
        any(grepl(rxn_check[i, Prod], info$EQUATION))
      if (!has_ko) return("noKOmatch")
      if (!has_cmpds) return("noCmpdmatch")
      "good"
    }
  })

  if (sum(rxn_id_check == "good") == 0) stop("No 'good' reaction entries. Check mapformula formatting.")
  rxn_check <- rxn_check[rxn_id_check == "good"]
  rxn_table <- merge(rxn_table, rxn_check, by = c("Rxn", "KO", "Reac", "Prod"), all.x = FALSE)

  # Compute stoichiometries per Reac/Prod from KEGG EQUATION
  goodrxns <- match(rxn_table[, Rxn], all_kegg$Reactions)
  stoichReac <- rep(1, nrow(rxn_table))
  stoichProd <- rep(1, nrow(rxn_table))

  for (j in seq_len(nrow(rxn_table))) {
    info <- all_kegg$Reaction_info[[goodrxns[j]]]
    if (!is.null(info$EQUATION)) {
      eqn <- strsplit(info$EQUATION, " <=> ", fixed = TRUE)
      comps <- unlist(lapply(eqn, strsplit, " \\+ "))
      coefs <- unique(comps[grepl(" ", comps)])
    } else {
      coefs <- character()
    }
    if (length(coefs) > 0) {
      sp <- strsplit(coefs, " ", fixed = TRUE)
      comp <- sapply(sp, function(x) x[2])
      reac_coef <- unlist(sp[comp == rxn_table[j, Reac]])[1]
      if (length(reac_coef) && !is.na(as.numeric(reac_coef))) stoichReac[j] <- as.numeric(reac_coef)
      prod_coef <- unlist(sp[comp == rxn_table[j, Prod]])[1]
      if (length(prod_coef) && !is.na(as.numeric(prod_coef))) stoichProd[j] <- as.numeric(prod_coef)
    }
  }

  rxn_table[, stoichReac := stoichReac]
  rxn_table[, stoichProd := stoichProd]
  rxn_table[, Path := gsub(" ", "", Path)]

  if (write_out) {
    stamp <- format(Sys.time(), "%Y-%m-%d_%H%M")
    fn <- paste0("full_network_template_all_info_", stamp, ".tsv")
    fwrite(rxn_table, fn, sep = "\t", quote = FALSE)
    message("Wrote network template -> ", fn)
  }
  rxn_table
}

