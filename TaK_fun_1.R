# 1. DINAMIC FUNCTION -----------------
calculate_TaK_shiny <- function(data, weight_vector) {
  meta_cols <- c("Dataset", "Abundance", "Row_Sum_W", "Row_Sum_P", "Total_N", 
                 "n_Obs", "n_Prop", "W_Obs", "W_Prop", "n_Lineages", "TR", "TC")
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
      n_Lineages = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      TR = ifelse(W_Prop > 0, W_Obs / W_Prop, 0),
      TC = ifelse(n_Prop > 0, n_Obs / n_Prop, 0)
    )
  return(main_results)
}