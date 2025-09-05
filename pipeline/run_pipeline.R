#!/usr/bin/env Rscript
# Run end-to-end: PICRUSt2 outputs -> PRMT scores -> basic QC plots

suppressPackageStartupMessages({
  source(file.path("R", "01_packages.R"))
  source(file.path("R", "02_io.R"))
  source(file.path("R", "03_preprocessing.R"))
  source(file.path("R", "04_prmt_mapping.R"))
  source(file.path("R", "05_statistics.R"))
  source(file.path("R", "06_visualization.R"))
  source(file.path("R", "99_utils.R"))
})

args <- commandArgs(trailingOnly = TRUE)
cfg_file <- if (length(args) >= 1) args[1] else file.path("pipeline", "config.yaml")
cfg <- yaml::read_yaml(cfg_file)

dir.create(cfg$output$dir, showWarnings = FALSE, recursive = TRUE)

# --- Load
abund <- read_picrust2_abundance(cfg$inputs$abundance_table, id_col = cfg$inputs$id_col)
meta  <- read_metadata(cfg$inputs$metadata)
maps  <- read_prmt_mapping(cfg$inputs$prmt_mapping)

# --- Align
aligned <- align_samples(abund, id_col = cfg$inputs$id_col, metadata = meta)
abund_a <- aligned$abund
meta_a  <- aligned$meta

# --- Normalize
if (identical(cfg$processing$normalize, "relative")) {
  abund_n <- rel_abundance(abund_a, id_col = cfg$inputs$id_col)
} else {
  abund_n <- abund_a
}

# --- Scores
prmt <- compute_prmt_scores(abund_n, id_col = cfg$inputs$id_col, mapping = maps, method = cfg$scoring$method)

# --- Optional z-score matrix for heatmaps
prmt_z <- zscore_prmt(prmt)

# --- Save outputs
out_scores <- file.path(cfg$output$dir, paste0(cfg$output$prefix, "_scores.tsv"))
out_z      <- file.path(cfg$output$dir, paste0(cfg$output$prefix, "_scores_z.tsv"))
readr::write_tsv(prmt, out_scores)
readr::write_tsv(prmt_z, out_z)

# --- Quick plots
p1 <- plot_prmt_density(prmt)
ggsave(file.path(cfg$output$dir, paste0(cfg$output$prefix, "_density.png")), p1, width = 7, height = 5, dpi = 150)

p2 <- plot_prmt_heatmap(prmt_z, top_n = 40)
ggsave(file.path(cfg$output$dir, paste0(cfg$output$prefix, "_heatmap.png")), p2, width = 8, height = 10, dpi = 150)

# --- Session info for reproducibility
writeLines(c(capture.output(sessionInfo())), file.path(cfg$output$dir, "sessionInfo.txt"))
message("Done. Results in: ", cfg$output$dir)

