# 1. DINAMIC FUNCTION -----------------
# NOTE on naming: `unique_taxa` counts distinct taxonomic lineage rows per
# Dataset x group_col (i.e. the "taxa uniques" count, analogous to the first
# number in each Table 2 cell). `Total_N` is the summed Abundance/count column
# (individuals or a density-style estimate depending on the source data), NOT
# a count of underlying occurrence records — the two are not interchangeable
# (see REPORT_revisao_R2.md, Task 1). This function does not compute a
# per-rank "terminal identification record" breakdown like manuscript
# Table 2; that is a structurally different, rank-by-rank calculation.
calculate_TaK_shiny <- function(data, weight_vector) {
  meta_cols <- c("Dataset", "Abundance", "Row_Sum_W", "Row_Sum_P", "Total_N",
                 "n_Obs", "n_Prop", "W_Obs", "W_Prop", "unique_taxa", "TR", "TC")
  tax_cols <- setdiff(names(data), meta_cols)
  
  if(length(tax_cols) == 0) return(NULL)
  
  weights <- if (length(weight_vector) != length(tax_cols)) 
    seq_along(tax_cols) else weight_vector
  rank_weights <- setNames(weights, tax_cols)
  max_potential_weight <- sum(weights)
  
  data$Abundance <- as.numeric(replace_na(as.character(data$Abundance), "0"))
  tax_matrix <- data[, tax_cols, drop = FALSE]
  is_identified <- !is.na(tax_matrix) & tax_matrix != "NA" & tax_matrix != ""
  
  data$Row_Sum_W <- as.matrix(is_identified) %*% rank_weights
  data$Row_Sum_P <- max_potential_weight
  
  group_col <- tax_cols[1]
  
  main_results <- data %>%
    group_by(Dataset, .data[[group_col]]) %>%
    summarise(
      n_Obs   = sum(Abundance * Row_Sum_W, na.rm = TRUE),
      n_Prop  = sum(Abundance * Row_Sum_P, na.rm = TRUE),
      W_Obs   = sum(Row_Sum_W, na.rm = TRUE),
      W_Prop  = sum(Row_Sum_P, na.rm = TRUE),
      Total_N = sum(Abundance, na.rm = TRUE),
      unique_taxa = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      TR = ifelse(W_Prop > 0, W_Obs / W_Prop, 0),
      TC = ifelse(n_Prop > 0, n_Obs / n_Prop, 0)
    )
  return(main_results)
}

# 2. SAMPLING-SENSITIVITY (RAREFACTION) FUNCTION -----------------
# Checks how stable TR and TC are as sampling effort (number of unique
# taxonomic lineages) decreases, via row-level bootstrap subsampling
# without replacement, computed per Dataset (not per group_col — this is a
# Dataset-wide TR/TC, matching how the manuscript's Table 2 regions are
# each a single Dataset). TR is an unweighted mean over lineages, so it is
# expected to stay ~flat across subsampling depth; TC is abundance-weighted
# and can swing sharply if a dataset's abundance is concentrated in a
# handful of coarsely-identified lineages (see
# REPORT_sampling_sensitivity.md for the empirical case that motivated
# this feature — one CCZ record alone held 72% of that region's abundance
# and moved TC 3-fold depending on whether it was sampled).
#
# progress_callback(step, total_steps), if given, is called once per
# (Dataset x fraction) combination — lets a caller drive a Shiny progress
# bar without this function depending on shiny itself.
rarefy_tak <- function(data, weight_vector, fractions = seq(0.1, 1, by = 0.1),
                        n_boot = 100, seed = 42, progress_callback = NULL) {
  meta_cols <- c("Dataset", "Abundance", "Row_Sum_W", "Row_Sum_P", "Total_N",
                 "n_Obs", "n_Prop", "W_Obs", "W_Prop", "unique_taxa", "TR", "TC")
  tax_cols <- setdiff(names(data), meta_cols)
  if (length(tax_cols) == 0) return(NULL)

  weights <- if (length(weight_vector) != length(tax_cols))
    seq_along(tax_cols) else weight_vector
  max_w <- sum(weights)

  data$Abundance <- as.numeric(replace_na(as.character(data$Abundance), "0"))
  tax_matrix <- data[, tax_cols, drop = FALSE]
  is_identified <- !is.na(tax_matrix) & tax_matrix != "NA" & tax_matrix != ""
  data$Row_Sum_W <- as.numeric(as.matrix(is_identified) %*% weights)

  subset_tr_tc <- function(d) {
    n <- nrow(d)
    if (n == 0 || sum(d$Abundance) == 0) return(c(TR = NA_real_, TC = NA_real_))
    TR <- sum(d$Row_Sum_W) / (n * max_w)
    TC <- sum(d$Abundance * d$Row_Sum_W) / (sum(d$Abundance) * max_w)
    c(TR = TR, TC = TC)
  }

  set.seed(seed)
  datasets <- unique(data$Dataset)
  total_steps <- length(datasets) * length(fractions)
  step <- 0
  results <- list()

  for (ds in datasets) {
    sub_full <- data[data$Dataset == ds, ]
    n_rows <- nrow(sub_full)
    if (n_rows < 4) next  # too few lineages for a meaningful rarefaction curve
    for (frac in fractions) {
      step <- step + 1
      if (!is.null(progress_callback)) progress_callback(step, total_steps)
      k <- min(n_rows, max(2, round(frac * n_rows)))
      boot <- t(sapply(seq_len(n_boot), function(i) {
        samp <- sub_full[sample.int(n_rows, k, replace = FALSE), ]
        subset_tr_tc(samp)
      }))
      results[[length(results) + 1]] <- data.frame(
        Dataset = ds, Fraction = frac, n_rows_total = n_rows, k = k,
        TR_mean = mean(boot[, "TR"], na.rm = TRUE), TR_sd = sd(boot[, "TR"], na.rm = TRUE),
        TC_mean = mean(boot[, "TC"], na.rm = TRUE), TC_sd = sd(boot[, "TC"], na.rm = TRUE)
      )
    }
  }
  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# Draws TR/TC-vs-subsampling-fraction curves (one line per Dataset), mean
# +/- 1 SD across bootstrap draws. Shared by the in-app plot and the
# PNG/SVG download handlers, so the three never diverge.
draw_rarefaction_plot <- function(rare_df) {
  plot_df <- rare_df %>%
    tidyr::pivot_longer(cols = c(TR_mean, TC_mean), names_to = "Metric", values_to = "Mean") %>%
    mutate(
      SD = ifelse(Metric == "TR_mean", TR_sd, TC_sd),
      Metric = recode(Metric, TR_mean = "Taxonomic Resolution (TR)",
                       TC_mean = "Taxonomic Completeness (TC)"),
      Metric = factor(Metric, levels = c("Taxonomic Resolution (TR)",
                                          "Taxonomic Completeness (TC)"))
    )

  ggplot(plot_df, aes(x = Fraction, y = Mean, color = Dataset, fill = Dataset)) +
    geom_ribbon(aes(ymin = pmax(0, Mean - SD), ymax = pmin(1, Mean + SD)),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6) +
    scale_x_continuous(labels = scales::percent, breaks = sort(unique(plot_df$Fraction))) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~Metric) +
    labs(x = "Subsampled proportion of unique lineages", y = "Index value") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}