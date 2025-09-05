#!/usr/bin/env Rscript
# End-to-end: fetch → KGML → mapformula → template → (optional) genomic matrices

suppressPackageStartupMessages({
  source("R/01_packages.R")
  source("R/07_kegg_fetch.R")
  source("R/11_kegg_info.R")
  source("R/08_network_template.R")
  source("R/09_genomic_network.R")
  source("R/10_prmt_cmp_scores.R")
  library(yaml)
  library(readr)
  library(data.table)
})

# --- config ---
cfg_path <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "pipeline/config.yaml"
cfg <- yaml::read_yaml(cfg_path)

out_dir <- if (!is.null(cfg$output$dir)) cfg$output$dir else "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

kgml_dir <- if (!is.null(cfg$kegg$kgml_dir)) cfg$kegg$kgml_dir else "kgml_files"
mapfile_nm <- if (!is.null(cfg$output$mapformula_file)) cfg$output$mapformula_file else "reaction_mapformula.tsv"
mapfile <- file.path(out_dir, mapfile_nm)
dir.create(dirname(mapfile), showWarnings = FALSE, recursive = TRUE)

# --- 1) Reactions / Pathways list ---
lists <- kegg_lists()

# --- 2) Optional reaction→pathway scan (resume-friendly) ---
if (isTRUE(cfg$kegg$scan_reactions)) {
  checkpoint <- if (!is.null(cfg$kegg$checkpoint)) cfg$kegg$checkpoint else "processed_reactions.RDS"
  unique_paths <- reaction_pathway_scan(lists$reaction_ids, checkpoint_file = checkpoint)
  saveRDS(unique_paths, file.path(out_dir, "unique_pathways.RDS"))
}

# --- 3) KGML download (skip if cache exists) ---
n_xml <- length(list.files(kgml_dir, pattern = "\\.xml$", full.names = TRUE))
if (n_xml > 0L) {
  message("KGML cache detected (", n_xml, " files) — skipping download.")
} else {
  download_kgml_ko(lists$pathway_ids, kgml_dir = kgml_dir, overwrite = FALSE)
}

# --- 4) Build reaction_mapformula from KGML ---
invisible(parse_kgml_to_mapformula(kgml_dir, out_file = mapfile))

# --- 5) Ensure KO↔Reaction map, load it ---
kos_rxns_path <- if (!is.null(cfg$kegg$kos_to_rxns_file)) cfg$kegg$kos_to_rxns_file else "inst/extdata/reference/kos_rxns.tsv"
ensure_kos_rxns(kos_rxns_path)  # download if missing

kos_dt <- data.table::fread(kos_rxns_path, header = FALSE, sep = "\t",
                            colClasses = "character", col.names = c("KO","Rxn"))

# --- 6) Build all_kegg s robustním REST fetch + checkpoint ---
kolist     <- if (!is.null(cfg$kegg$kolist)) cfg$kegg$kolist else character()
all_pref   <- if (!is.null(cfg$kegg$all_kegg_prefix)) cfg$kegg$all_kegg_prefix else "all_kegg"
out_prefix <- file.path(out_dir, all_pref)
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

checkpoint_file <- file.path(out_dir, paste0(all_pref, "_fetch_checkpoint.rds"))

all_kegg <- get_kegg_reaction_info(
  kos_to_rxns_file = kos_dt,
  kolist           = kolist,
  save_out         = TRUE,
  out_prefix       = out_prefix,
  checkpoint_file  = checkpoint_file,
  chunk_size       = 5,     # malé dávky pro REST
  sleep_sec        = 1,
  long_sleep_sec   = 60,
  max_retries      = 6
)

# --- 7) Community network template (Rxn, KO, Reac, Prod, stoich) ---
template <- generate_network_template_kegg(
  mapformula_file = mapfile,
  all_kegg        = all_kegg,
  write_out       = TRUE
)

message("PRMT mapping build: DONE.")

