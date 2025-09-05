# Generate stoichiometric and adjacency-like matrices per KO subset
# Refactor of your generate_genomic_network()

suppressPackageStartupMessages({
  library(data.table)
})

#' Generate genomic network for a set of KOs using KEGG template or MetaCyc
#' @param kos character KOs to include
#' @param keggSource "KeggTemplate" | "loadNet" | "metacyc"
#' @param degree_filter numeric; remove currency metabolites above this degree (0 = keep all)
#' @param minpath_file optional TSV with column 'Path' to restrict to minimal pathways
#' @param normalize logical; row-normalize negative and positive parts separately
#' @param rxn_table KEGG template (output of generate_network_template_kegg)
#' @param networkFile RData path when keggSource == "loadNet"
#' @return list(network_mat, stoich_mat, rxn_table_subset)
generate_genomic_network <- function(kos,
                                     keggSource = "KeggTemplate",
                                     degree_filter = 0,
                                     minpath_file = "",
                                     normalize = TRUE,
                                     rxn_table = NULL,
                                     networkFile = "") {
  if (!keggSource %in% c("loadNet", "KeggTemplate", "metacyc"))
    stop("Invalid keggSource. Use 'loadNet', 'KeggTemplate', or 'metacyc'.")

  if (keggSource == "loadNet") {
    load(networkFile) # should define 'allnet'
    return(allnet)
  }

  if (keggSource == "KeggTemplate") {
    if (is.null(rxn_table)) stop("Provide rxn_table (KEGG template).")

    tbl <- rxn_table[KO %in% kos]
    if (nrow(tbl) == 0) stop("No reactions found for the given KO set.")

    if (nzchar(minpath_file)) {
      minpaths <- fread(minpath_file, colClasses = "character")
      setnames(minpaths, "Path")
      tbl2 <- tbl[Path %in% minpaths[, Path]]
      tbl_extra <- tbl[!(Rxn %in% tbl2[, Rxn])]
      tbl <- rbind(tbl2, tbl_extra)
    }

    # Drop unused columns
    for (dropcol in c("Path", "ReacProd", "Rxn")) {
      if (dropcol %in% names(tbl)) tbl[, (dropcol) := NULL]
    }
    setkey(tbl, NULL)
    tbl <- unique(tbl)

    cmpds <- unique(c(tbl[, Prod], tbl[, Reac]))

    if (degree_filter != 0) {
      degree <- sapply(cmpds, function(x) tbl[Prod == x | Reac == x, length(unique(KO))])
      keep <- degree < degree_filter
      cmpds <- cmpds[keep]
      tbl  <- tbl[(Prod %in% cmpds) & (Reac %in% cmpds)]
      cmpds <- unique(c(tbl[, Prod], tbl[, Reac]))
    }

    goodkos <- unique(tbl[, KO])

    # Build stoichiometric matrices (compounds x KOs)
    network_mat <- matrix(0, nrow = length(cmpds), ncol = length(goodkos))
    stoich_mat  <- matrix(NA_real_, nrow = length(cmpds), ncol = length(goodkos))

    for (j in seq_len(nrow(tbl))) {
      i_reac <- match(tbl[j, Reac], cmpds)
      i_prod <- match(tbl[j, Prod], cmpds)
      i_ko   <- match(tbl[j, KO],   goodkos)

      network_mat[i_reac, i_ko] <- network_mat[i_reac, i_ko] - tbl[j, stoichReac]
      network_mat[i_prod, i_ko] <- network_mat[i_prod, i_ko] + tbl[j, stoichProd]

      stoich_mat[i_reac, i_ko] <- if (is.na(stoich_mat[i_reac, i_ko])) -tbl[j, stoichReac] else stoich_mat[i_reac, i_ko] - tbl[j, stoichReac]
      stoich_mat[i_prod, i_ko] <- if (is.na(stoich_mat[i_prod, i_ko]))  tbl[j, stoichProd] else stoich_mat[i_prod, i_ko] + tbl[j, stoichProd]
    }

    if (normalize) {
      negsums <- apply(network_mat, 1, function(x) sum(x[x < 0]))
      possums <- apply(network_mat, 1, function(x) sum(x[x > 0]))
      for (r in seq_len(nrow(network_mat))) {
        negkos <- which(network_mat[r, ] < 0)
        poskos <- which(network_mat[r, ] > 0)
        if (length(negkos)) network_mat[r, negkos] <- network_mat[r, negkos] / abs(negsums[r])
        if (length(poskos)) network_mat[r, poskos] <- network_mat[r, poskos] /  possums[r]
      }
    }

    network_mat <- as.data.frame(network_mat, stringsAsFactors = FALSE)
    colnames(network_mat) <- goodkos
    rownames(network_mat) <- cmpds

    stoich_mat <- as.data.frame(stoich_mat, stringsAsFactors = FALSE)
    colnames(stoich_mat) <- goodkos
    rownames(stoich_mat) <- cmpds

    return(list(network_mat, stoich_mat, tbl))
  }

  if (keggSource == "metacyc") {
    stop("MetaCyc branch preserved in original script; port if/when needed.")
  }
}

