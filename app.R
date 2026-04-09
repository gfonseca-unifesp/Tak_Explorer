# --- 0. INSTALAÇÃO AUTOMATIZADA DE PACOTES ---
install_missing_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    message(paste("Instalando pacotes ausentes:", paste(new_packages, collapse = ", ")))
    install.packages(new_packages, dependencies = TRUE)
  }
}

# Lista de pacotes necessários para o app
required_libs <- c("shiny", "bslib", "dplyr", "tidyr", "ggplot2", 
                   "ggrepel", "DT", "scales", "readr")

# Executa a verificação/instalação
install_missing_packages(required_libs)

# Carrega as bibliotecas
library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(DT)
library(scales)
library(readr)

# --- 1. FUNÇÃO DE CÁLCULO DINÂMICA ---
calculate_TaK_shiny <- function(data, weight_vector) {
  meta_cols <- c("Dataset", "Abundance", "Row_Sum_W", "Row_Sum_P", "Total_N", 
                 "n_Obs", "n_Prop", "W_Obs", "W_Prop", "n_Lineages", "TR", "TC")
  tax_cols <- setdiff(names(data), meta_cols)
  
  if(length(tax_cols) == 0) return(NULL)
  
  weights <- if (length(weight_vector) != length(tax_cols)) seq_along(tax_cols) else weight_vector
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

# --- 2. INTERFACE (UI) ---
ui <- page_navbar(
  title = "TaK Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("Data Editor",
            layout_sidebar(
              sidebar = sidebar(
                title = "Settings",
                fileInput("upload_csv", "Upload CSV", accept = c(".csv", ".txt")),
                hr(),
                textInput("weight_vec_str", "Weight Vector:", value = "1, 2, 3, 4, 5, 6"),
                numericInput("row_count", "Number of Lineages:", value = 10, min = 1),
                actionButton("add_column", "Add taxonomic column", icon = icon("plus")),
                hr(),
                selectInput("label_choice", "Point label (Biplot):", 
                            choices = c("Dataset" = "Dataset", "Lineage" = "Lineage")),
                actionButton("reset_data", "Reset App", class = "btn-warning w-100")
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  "Lineage Editor",
                  downloadButton("download_data", "Export CSV", class = "btn-sm")
                ),
                DTOutput("editable_table")
              )
            )
  ),
  
  nav_panel("Visualization",
            layout_column_wrap(
              width = 1/2,
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  "Biplot (TR vs TC)",
                  downloadButton("download_biplot", "PNG", class = "btn-sm")
                ),
                plotOutput("biplot", height = "500px")
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  "Frequency per Quadrant",
                  downloadButton("download_quadrant", "PNG", class = "btn-sm")
                ),
                plotOutput("quadrant_plot", height = "500px")
              )
            )
  ),
  
  nav_panel("Summary",
            card(
              card_header(
                class = "d-flex justify-content-between align-items-center",
                "Dataset Summary (Aggregated Indices)",
                downloadButton("download_summary", "Download Summary CSV", class = "btn-success btn-sm")
              ),
              DTOutput("summary_table")
            )
  )
)

# --- 3. SERVIDOR (SERVER) ---
server <- function(input, output, session) {
  
  v <- reactiveValues(data = data.frame(
    Dataset = rep("Exemplo", 5),
    Phylum = "Arthropoda", Class = "Insecta", Order = "Diptera", 
    Family = "Muscidae", Genus = "Musca", Species = "domestica",
    Abundance = c(10, 20, 30, 40, 50),
    stringsAsFactors = FALSE
  ))
  
  observeEvent(input$add_column, {
    new_col_name <- paste0("NewRank_", ncol(v$data) - 1)
    abund_idx <- which(names(v$data) == "Abundance")
    v$data <- v$data %>% mutate(!!new_col_name := "") %>% select(1:(abund_idx-1), last_col(), Abundance)
    updateTextInput(session, "weight_vec_str", value = paste0(input$weight_vec_str, ", ", ncol(v$data)-2))
  })
  
  current_weights <- reactive({ as.numeric(unlist(strsplit(input$weight_vec_str, ","))) })
  
  observeEvent(input$upload_csv, {
    req(input$upload_csv)
    line <- readLines(input$upload_csv$datapath, n = 1)
    sep <- if(grepl(";", line)) ";" else ","
    df <- read_delim(input$upload_csv$datapath, delim = sep, show_col_types = FALSE)
    
    found_n <- intersect(names(df), c("individualCount", "n", "Abundance", "abundance", "count"))[1]
    if(!is.na(found_n)) df <- df %>% rename(Abundance = !!sym(found_n))
    if(!"Abundance" %in% names(df)) df$Abundance <- 1
    if(!"Dataset" %in% names(df)) df$Dataset <- "Uploaded"
    
    tax_cols <- setdiff(names(df), c("Dataset", "Abundance"))
    v$data <- df %>% select(Dataset, all_of(tax_cols), Abundance) %>% as.data.frame()
    updateNumericInput(session, "row_count", value = nrow(v$data))
    updateTextInput(session, "weight_vec_str", value = paste(seq_along(tax_cols), collapse = ", "))
  })
  
  observeEvent(input$row_count, {
    req(v$data)
    if (input$row_count > nrow(v$data)) {
      extra <- v$data[rep(1, input$row_count - nrow(v$data)), ]
      extra[extra != ""] <- ""; extra$Abundance <- 1; extra$Dataset <- "Novo"
      v$data <- rbind(v$data, extra)
    } else { v$data <- v$data[1:input$row_count, ] }
  })
  
  output$editable_table <- renderDT({
    datatable(v$data, editable = list(target = 'all'), options = list(pageLength = 10, dom = 'tip'))
  })
  
  observeEvent(input$editable_table_cell_edit, { v$data <<- editData(v$data, input$editable_table_cell_edit) })
  
  processed_results <- reactive({ req(v$data); calculate_TaK_shiny(v$data, current_weights()) })
  
  summary_df <- reactive({
    res <- processed_results()
    req(res)
    res %>%
      group_by(Dataset) %>%
      summarise(
        Mean_TR = mean(TR, na.rm = TRUE),
        Mean_TC = mean(TC, na.rm = TRUE),
        Total_Individuals = sum(Total_N, na.rm = TRUE),
        Taxa_Groups_Count = n(),
        Total_Lineages = sum(n_Lineages, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(across(c(Mean_TR, Mean_TC), ~round(., 3)))
  })
  
  output$summary_table <- renderDT({
    datatable(summary_df(), options = list(pageLength = 10, dom = 't'))
  })
  
  plot_biplot <- reactive({
    res <- processed_results(); req(res)
    label_col <- if(input$label_choice == "Dataset") names(res)[1] else names(res)[2]
    ggplot(res, aes(x = TC, y = TR, color = Dataset)) +
      annotate("rect", xmin=0.5, xmax=1, ymin=0.5, ymax=1, fill="#2ecc71", alpha=0.1) +
      annotate("rect", xmin=0, xmax=0.5, ymin=0.5, ymax=1, fill="#3498db", alpha=0.1) +
      annotate("rect", xmin=0, xmax=0.5, ymin=0, ymax=0.5, fill="#e74c3c", alpha=0.1) +
      annotate("rect", xmin=0.5, xmax=1, ymin=0, ymax=0.5, fill="#f1c40f", alpha=0.1) +
      geom_point(aes(size=Total_N), alpha=0.7) +
      geom_text_repel(aes(label = .data[[label_col]])) +
      geom_abline(intercept=0, slope=1, linetype="dashed") +
      scale_x_continuous(limits=c(0, 1.05)) + scale_y_continuous(limits=c(0, 1.05)) +
      theme_minimal() + theme(legend.position = "bottom")
  })
  
  plot_quadrant <- reactive({
    res <- processed_results(); req(res)
    quad_data <- res %>%
      mutate(Quadrant = case_when(
        TC >= 0.5 & TR >= 0.5 ~ "Well resolved",
        TC <  0.5 & TR >= 0.5 ~ "Under resolved",
        TC <  0.5 & TR <  0.5 ~ "Deficient",
        TC >= 0.5 & TR <  0.5 ~ "Abundance biased"
      )) %>%
      group_by(Dataset, Quadrant) %>%
      summarise(Count = n(), .groups = 'drop_last') %>%
      mutate(Perc = Count / sum(Count))
    
    ggplot(quad_data, aes(x = Dataset, y = Perc, fill = Quadrant)) +
      geom_bar(stat = "identity", position = "fill", color = "white") +
      geom_text(aes(label = percent(Perc, accuracy = 1)), position = position_fill(vjust = 0.5), fontface="bold") +
      scale_fill_manual(values = c("Well resolved"="#2ecc71", "Under resolved"="#3498db", "Abundance biased"="#f1c40f", "Deficient"="#e74c3c")) +
      scale_y_continuous(labels = percent) +
      theme_minimal() + theme(legend.position = "bottom")
  })
  
  output$biplot <- renderPlot({ plot_biplot() })
  output$quadrant_plot <- renderPlot({ plot_quadrant() })
  
  output$download_data <- downloadHandler(
    filename = function() { paste("tak_raw_data_", Sys.Date(), ".csv", sep="") },
    content = function(file) { write.csv(v$data, file, row.names = FALSE) }
  )
  
  output$download_summary <- downloadHandler(
    filename = function() { paste("tak_summary_", Sys.Date(), ".csv", sep="") },
    content = function(file) { write.csv(summary_df(), file, row.names = FALSE) }
  )
  
  output$download_biplot <- downloadHandler(
    filename = function() { paste("biplot_", Sys.Date(), ".png", sep="") },
    content = function(file) { ggsave(file, plot = plot_biplot(), width = 10, height = 7, dpi = 300) }
  )
  
  output$download_quadrant <- downloadHandler(
    filename = function() { paste("quadrants_", Sys.Date(), ".png", sep="") },
    content = function(file) { ggsave(file, plot = plot_quadrant(), width = 10, height = 7, dpi = 300) }
  )
  
  observeEvent(input$reset_data, { session$reload() })
}

shinyApp(ui, server)