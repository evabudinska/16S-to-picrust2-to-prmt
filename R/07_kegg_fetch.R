# 07_kegg_fetch.R — KEGG fetching & KGML parsing to reaction_mapformula

suppressPackageStartupMessages({
  library(KEGGREST)
  library(pathview)   # for download.kegg()
  library(XML)
  library(data.table)
})

#' Get KEGG reactions and KO pathways (robust id parsing)
#' @return list(reaction_ids, pathway_ids, pathways_df)
kegg_lists <- function() {
  reaction_list <- KEGGREST::keggList("reaction")
  reaction_ids  <- names(reaction_list)

  pathway_list <- KEGGREST::keggList("pathway/ko")
  raw_ids <- names(pathway_list)
  # extract last 5 digits (e.g., "path:ko01100" -> "01100", "ko01100" -> "01100")
  five <- regmatches(raw_ids, regexpr("[0-9]{5}$", raw_ids))
  pathway_ids <- unique(five[nzchar(five)])

  if (!length(pathway_ids)) stop("No KO pathway IDs parsed from keggList('pathway/ko').")

  pathways_df <- data.frame(
    ID   = pathway_ids,
    Name = as.character(pathway_list)[seq_along(pathway_list)],
    stringsAsFactors = FALSE
  )

  list(
    reaction_ids = reaction_ids,
    pathway_ids  = pathway_ids,
    pathways_df  = pathways_df
  )
}

#' Resume-friendly traversal of reaction -> pathway links via KEGGREST::keggGet
reaction_pathway_scan <- function(reaction_ids,
                                  checkpoint_file = "processed_reactions.RDS",
                                  sleep_sec = 1) {
  state <- if (file.exists(checkpoint_file)) {
    readRDS(checkpoint_file)
  } else {
    list(index = 1L, all_pathways = character())
  }

  for (i in seq.int(state$index, length(reaction_ids))) {
    rid <- reaction_ids[i]
    try({
      rxn <- KEGGREST::keggGet(rid)[[1]]
      if (!is.null(rxn$PATHWAY)) {
        rn_ids <- names(rxn$PATHWAY)
        rn_ids <- gsub("^rn", "", rn_ids)
        state$all_pathways <- c(state$all_pathways, rn_ids)
      }
      state$index <- i + 1L
      saveRDS(state, checkpoint_file)
    }, silent = TRUE)
    Sys.sleep(sleep_sec)
  }
  unique(state$all_pathways)
}

#' Download KGML files for KO pathways
#' Accepts "01100" or "ko01100", saves as ko01100.xml
download_kgml_ko <- function(pathway_ids, kgml_dir = "kgml_files", overwrite = FALSE) {
  dir.create(kgml_dir, showWarnings = FALSE, recursive = TRUE)

  ids <- gsub("^path:ko", "", pathway_ids)
  ids <- gsub("^ko", "", ids)
  ids <- ids[grepl("^[0-9]{5}$", ids)]
  if (!length(ids)) return(character())

  existing <- list.files(kgml_dir, pattern = "^ko[0-9]{5}\\.xml$", full.names = FALSE)
  have_ids <- sub("^ko|\\.xml$", "", existing)

  todo_ids <- if (isTRUE(overwrite)) ids else setdiff(ids, have_ids)
  if (!length(todo_ids)) return(sort(unique(have_ids)))

  message("Downloading ", length(todo_ids), " KGML files to '", kgml_dir, "' ...")
  tryCatch({
    pathview::download.kegg(
      pathway.id = todo_ids,
      species    = "ko",
      file.type  = "xml",
      kegg.dir   = kgml_dir
    )
  }, error = function(e) {
    stop("KGML download failed: ", conditionMessage(e))
  })

  sort(unique(c(have_ids, todo_ids)))
}

#' Parse KGML files to reaction_mapformula lines
#' Output lines: "Rxxxxx: ppppp: C00001 + C00002 <=> C00003"
parse_kgml_to_mapformula <- function(kgml_dir = "kgml_files", out_file = NULL) {
  files <- list.files(kgml_dir, pattern = "\\.xml$", full.names = TRUE)
  if (!length(files)) stop("No KGML files found in: ", kgml_dir)

  lines_out <- character()

  for (file in files) {
    doc <- XML::xmlParse(file)
    root <- XML::xmlRoot(doc)

    fname <- basename(file)
    pathway_id <- sub("^ko", "", sub("\\.xml$", "", fname))  # "ko01100.xml" -> "01100"

    reactions  <- XML::getNodeSet(root, "//reaction")
    if (!length(reactions)) next

    for (rxn in reactions) {
      reaction_names <- XML::xmlGetAttr(rxn, "name")
      reaction_ids   <- gsub("^(rn:|gl:)", "", unlist(strsplit(reaction_names, " ")))

      direction  <- XML::xmlGetAttr(rxn, "type", default = "reversible")
      dir_symbol <- ifelse(direction == "reversible", "<=>", "=>")

      subs  <- XML::xpathSApply(rxn, "substrate", XML::xmlGetAttr, "name")
      prods <- XML::xpathSApply(rxn, "product",  XML::xmlGetAttr, "name")
      if (!length(subs) || !length(prods)) next

      subs_str  <- paste(sub("^cpd:", "", subs),  collapse = " + ")
      prods_str <- paste(sub("^cpd:", "", prods), collapse = " + ")

      for (rid in reaction_ids) {
        lines_out <- c(
          lines_out,
          paste0(rid, ": ", pathway_id, ": ", subs_str, " ", dir_symbol, " ", prods_str)
        )
      }
    }
  }

  lines_out <- sort(unique(lines_out))
  if (!is.null(out_file)) {
    writeLines(lines_out, out_file)
    message("Wrote reaction_mapformula with ", length(lines_out), " entries -> ", out_file)
  }
  lines_out
}

#' Convenience: end-to-end fetch+parse
kegg_fetch_and_parse <- function(kgml_dir = "kgml_files",
                                 out_file = NULL,
                                 subset_ids = NULL) {
  ids <- if (is.null(subset_ids)) kegg_lists()$pathway_ids else subset_ids
  download_kgml_ko(ids, kgml_dir = kgml_dir, overwrite = FALSE)
  parse_kgml_to_mapformula(kgml_dir = kgml_dir, out_file = out_file)
}

