#' Minimal visualizations for quick QC

#' Density plot of PRMT scores for a sample (or overlay)
plot_prmt_density <- function(prmt_wide) {
  long <- prmt_wide |>
    tidyr::pivot_longer(-prmt_id, names_to = "sample_id", values_to = "score")
  ggplot(long, aes(score, group = sample_id)) +
    geom_density(alpha = 0.2) +
    labs(x = "PRMT score", y = "Density") +
    theme_bw()
}

#' Simple heatmap using ggplot (no ComplexHeatmap dependency)
plot_prmt_heatmap <- function(prmt_wide, top_n = 50) {
  # Select top variable PRMTs for readability
  scols <- setdiff(names(prmt_wide), "prmt_id")
  vars  <- prmt_wide |>
    mutate(var = apply(across(all_of(scols)), 1, var, na.rm = TRUE)) |>
    arrange(desc(var)) |>
    slice(seq_len(min(top_n, n()))) |>
    select(-var)

  long <- vars |>
    tidyr::pivot_longer(-prmt_id, names_to = "sample_id", values_to = "score")

  ggplot(long, aes(sample_id, prmt_id, fill = score)) +
    geom_tile() +
    scale_x_discrete(position = "top") +
    labs(x = "Sample", y = "PRMT", fill = "Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0))
}

