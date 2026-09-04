# 1. CORE FUNCTION -----------------
# NOTE on naming: `unique_taxa` counts distinct taxonomic lineage rows per
# Group x Rank (i.e. the "taxa uniques" count, analogous to the first
# number in each Table 2 cell). `Total_N` is the summed Abundance/count column
# (individuals or a density-style estimate depending on the source data), NOT
# a count of underlying occurrence records — the two are not interchangeable
# (see REPORT_revisao_R2.md, Task 1). This function does not compute a
# per-rank "terminal identification record" breakdown like manuscript
# Table 2; that is a structurally different, rank-by-rank calculation.
#
# GENERALIZATION: this used to always group by (Dataset, first taxonomic
# rank column). It now groups by (primary_col, rank_col), both selectable:
# - primary_col: the "compare across" dimension shown in Summary/Biplot/
#   quadrant plot/Sampling Sensitivity -- "Dataset" by default, or "Year"
#   when that optional column is present in the data.
# - rank_col: which taxonomic rank column to aggregate by -- defaults to
#   the first one (e.g. Phylum, matching the original behavior exactly)
#   but can be any rank column, which is what powers the "Rank Comparison"
#   tab (grouping by Genus, Species, etc. instead).
# The two grouping columns are always renamed to `Group` and `Rank` in the
# output, regardless of which real column fed them, so every downstream
# consumer (plots, tables, the rarefaction tab) has one stable pair of
# names to depend on instead of guessing column positions/names.
calculate_TaK_shiny <- function(data, weight_vector, primary_col = "Dataset", rank_col = NULL) {
  meta_cols <- c("Dataset", "Year", "Abundance", "Row_Sum_W", "Row_Sum_P", "Total_N",
                 "n_Obs", "n_Prop", "W_Obs", "W_Prop", "unique_taxa", "TR", "TC")
  tax_cols <- setdiff(names(data), meta_cols)

  if (length(tax_cols) == 0) return(NULL)
  if (is.null(rank_col) || !(rank_col %in% tax_cols)) rank_col <- tax_cols[1]
  if (!(primary_col %in% names(data))) primary_col <- "Dataset"

  weights <- if (length(weight_vector) != length(tax_cols))
    seq_along(tax_cols) else weight_vector
  rank_weights <- setNames(weights, tax_cols)
  max_potential_weight <- sum(weights)

  data$Abundance <- as.numeric(replace_na(as.character(data$Abundance), "0"))
  tax_matrix <- data[, tax_cols, drop = FALSE]
  is_identified <- !is.na(tax_matrix) & tax_matrix != "NA" & tax_matrix != ""

  data$Row_Sum_W <- as.matrix(is_identified) %*% rank_weights
  data$Row_Sum_P <- max_potential_weight

  main_results <- data %>%
    group_by(.data[[primary_col]], .data[[rank_col]]) %>%
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
  names(main_results)[1:2] <- c("Group", "Rank")
  # Group (Dataset or Year) and Rank (a taxon name) must always plot as
  # discrete categories, never as a continuous numeric axis/gradient --
  # Year in particular would otherwise come through as an integer and get
  # treated as continuous (fractional axis ticks, a color gradient instead
  # of a discrete legend) by every ggplot that maps it.
  main_results$Group <- as.character(main_results$Group)
  main_results$Rank <- as.character(main_results$Rank)
  return(main_results)
}

# 2. SAMPLING-SENSITIVITY (RAREFACTION) FUNCTION -----------------
# Checks how stable TR and TC are as sampling effort (number of unique
# taxonomic lineages) decreases, via row-level bootstrap subsampling
# without replacement, computed per primary_col value (Dataset or Year --
# same "compare across" dimension used everywhere else in the app; a
# Group-wide TR/TC, matching how the manuscript's Table 2 regions are
# each a single Dataset). TR is an unweighted mean over lineages, so it is
# expected to stay ~flat across subsampling depth; TC is abundance-weighted
# and can swing sharply if a Group's abundance is concentrated in a
# handful of coarsely-identified lineages (see
# REPORT_sampling_sensitivity.md for the empirical case that motivated
# this feature — one CCZ record alone held 72% of that region's abundance
# and moved TC 3-fold depending on whether it was sampled).
#
# progress_callback(step, total_steps), if given, is called once per
# (Group x fraction) combination — lets a caller drive a Shiny progress
# bar without this function depending on shiny itself.
rarefy_tak <- function(data, weight_vector, primary_col = "Dataset",
                        fractions = seq(0.1, 1, by = 0.1),
                        n_boot = 100, seed = 42, progress_callback = NULL) {
  meta_cols <- c("Dataset", "Year", "Abundance", "Row_Sum_W", "Row_Sum_P", "Total_N",
                 "n_Obs", "n_Prop", "W_Obs", "W_Prop", "unique_taxa", "TR", "TC")
  tax_cols <- setdiff(names(data), meta_cols)
  if (length(tax_cols) == 0) return(NULL)
  if (!(primary_col %in% names(data))) primary_col <- "Dataset"

  weights <- if (length(weight_vector) != length(tax_cols))
    seq_along(tax_cols) else weight_vector
  max_w <- sum(weights)

  data$Abundance <- as.numeric(replace_na(as.character(data$Abundance), "0"))
  tax_matrix <- data[, tax_cols, drop = FALSE]
  is_identified <- !is.na(tax_matrix) & tax_matrix != "NA" & tax_matrix != ""
  data$Row_Sum_W <- as.numeric(as.matrix(is_identified) %*% weights)
  data$Group <- as.character(data[[primary_col]])

  subset_tr_tc <- function(d) {
    n <- nrow(d)
    if (n == 0 || sum(d$Abundance) == 0) return(c(TR = NA_real_, TC = NA_real_))
    TR <- sum(d$Row_Sum_W) / (n * max_w)
    TC <- sum(d$Abundance * d$Row_Sum_W) / (sum(d$Abundance) * max_w)
    c(TR = TR, TC = TC)
  }

  set.seed(seed)
  groups <- unique(data$Group)
  total_steps <- length(groups) * length(fractions)
  step <- 0
  results <- list()

  for (g in groups) {
    sub_full <- data[data$Group == g, ]
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
        Group = g, Fraction = frac, n_rows_total = n_rows, k = k,
        TR_mean = mean(boot[, "TR"], na.rm = TRUE), TR_sd = sd(boot[, "TR"], na.rm = TRUE),
        TC_mean = mean(boot[, "TC"], na.rm = TRUE), TC_sd = sd(boot[, "TC"], na.rm = TRUE)
      )
    }
  }
  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# Draws TR/TC-vs-subsampling-fraction curves (one line per Group), mean
# +/- 1 SD across bootstrap draws. Shared by the in-app plot and the
# PNG/SVG download handlers, so the three never diverge. `group_label`
# only affects the legend title (e.g. "Dataset" or "Year") -- the data
# column is always called Group.
draw_rarefaction_plot <- function(rare_df, group_label = "Group") {
  plot_df <- rare_df %>%
    tidyr::pivot_longer(cols = c(TR_mean, TC_mean), names_to = "Metric", values_to = "Mean") %>%
    mutate(
      SD = ifelse(Metric == "TR_mean", TR_sd, TC_sd),
      Metric = recode(Metric, TR_mean = "Taxonomic Resolution (TR)",
                       TC_mean = "Taxonomic Completeness (TC)"),
      Metric = factor(Metric, levels = c("Taxonomic Resolution (TR)",
                                          "Taxonomic Completeness (TC)"))
    )

  ggplot(plot_df, aes(x = Fraction, y = Mean, color = Group, fill = Group)) +
    geom_ribbon(aes(ymin = pmax(0, Mean - SD), ymax = pmin(1, Mean + SD)),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6) +
    scale_x_continuous(labels = scales::percent, breaks = sort(unique(plot_df$Fraction))) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~Metric) +
    labs(x = "Subsampled proportion of unique lineages", y = "Index value",
         color = group_label, fill = group_label) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

# 3. RANK-COMPARISON BUBBLE PLOT -----------------
# Draws TR or TC (one metric at a time, mapped to both bubble size and
# color) with the chosen taxonomic rank's values on the Y axis and the
# app's current "compare by" dimension (Dataset or Year) on the X axis.
# `rank_df` is the output of calculate_TaK_shiny() called with rank_col
# set to whatever rank the user picked -- so a row here is one (Group,
# Rank-value) pair, e.g. (CCZ, Copepoda) with its own TR/TC.
#
# top_n limits the number of distinct Rank values shown (ranked by
# Total_N summed across all Group values, so the same set of rows is kept
# for every column -- a ragged Y axis between columns would be far harder
# to read). NULL/Inf shows everything.
draw_rank_comparison_plot <- function(rank_df, metric = "TR", group_label = "Group",
                                       rank_label = "Rank", top_n = 30) {
  if (is.finite(top_n) && top_n > 0) {
    keep <- rank_df %>%
      group_by(Rank) %>%
      summarise(total = sum(Total_N, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = top_n) %>%
      pull(Rank)
    rank_df <- rank_df %>% filter(Rank %in% keep)
  }

  metric_label <- if (metric == "TC") "Taxonomic Completeness (TC)" else "Taxonomic Resolution (TR)"

  ggplot(rank_df, aes(x = Group, y = Rank, size = .data[[metric]], color = .data[[metric]])) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(1, 10), limits = c(0, 1), name = metric_label) +
    scale_color_viridis_c(limits = c(0, 1), name = metric_label) +
    labs(x = group_label, y = rank_label,
         title = paste0(metric_label, " by ", rank_label, ", across ", group_label)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA))
}

# How many distinct Rank values a call to draw_rank_comparison_plot() with
# this top_n will actually end up showing -- shared by the screen render
# (dynamic plotOutput height) and the PNG download (dynamic ggsave height)
# so a deep rank with many categories (Genus, Species) always gets enough
# vertical room per row instead of a fixed plot height that made rows
# bleed into each other once there were more than ~15-20 of them (caught
# by rendering the real ISA_DeepData_2026.csv with Genus and looking at
# the actual PNG, not assumed from the code).
rank_comparison_n_shown <- function(rank_df, top_n = 30) {
  n_total <- length(unique(rank_df$Rank))
  if (is.finite(top_n) && top_n > 0) min(top_n, n_total) else n_total
}

rank_comparison_height_px <- function(n_shown) max(420, 140 + n_shown * 24)
rank_comparison_height_in <- function(n_shown) max(4.5, 1.4 + n_shown * 0.24)
