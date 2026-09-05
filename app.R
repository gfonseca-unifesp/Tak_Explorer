# --- 0. Install ---
install_missing_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    message(paste("Install absent packages:", paste(new_packages, collapse = ", ")))
    install.packages(new_packages, dependencies = TRUE)
  }
}

required_libs <- c("shiny", "bslib", "dplyr", "tidyr", "ggplot2", 
                   "ggrepel", "DT", "scales", "readr")

install_missing_packages(required_libs)

library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(DT)
library(scales)
library(readr)

# --- 1. Functions ---
# calculate_TaK_shiny() (main index), rarefy_tak() and draw_rarefaction_plot()
# (Sampling Sensitivity tab) all live in TaK_fun_1.R now, sourced from here
# instead of duplicated inline -- keeps the Shiny app and the standalone R
# function permanently in sync (see README, "Repository structure").
source("TaK_fun_1.R")

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
                numericInput("row_count", "Number of Lineages:", value = 66, min = 1),
                actionButton("add_column", "Add taxonomic column", icon = icon("plus")),
                actionButton("add_year_column", "Add Year column", icon = icon("plus")),
                hr(),
                selectInput("primary_col", "Compare by:", choices = c("Dataset" = "Dataset")),
                helpText("Governs Summary, Sampling Sensitivity and the ",
                         "Rank Comparison tab's X axis, plus whichever of ",
                         "\"Color points by\" (Visualization tab) is set ",
                         "to Dataset/Year. \"Year\" only appears once a ",
                         "Year column exists with at least one value ",
                         "filled in -- and when selected, records are ",
                         "pooled across all Datasets for each year (not ",
                         "cross-tabulated with Dataset)."),
                hr(),
                numericInput("low_confidence_threshold",
                             "Low-confidence threshold (min. records per group):",
                             value = 10, min = 0),
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
            layout_sidebar(
              sidebar = sidebar(
                title = "View settings",
                selectInput("viz_color_by", "Color points by:",
                            choices = c("Dataset / Year (Compare by)" = "Group",
                                        "Taxonomic level" = "Rank")),
                selectInput("viz_rank_col", "Taxonomic level:", choices = NULL),
                helpText("Which taxonomic rank column each point ",
                         "represents -- e.g. Phylum for the ISA data, but ",
                         "any rank your data has. Only matters for the ",
                         "two plots below (Summary always uses the first ",
                         "taxonomic column, unaffected by this)."),
                hr(),
                selectInput("viz_group_filter", "Show:", choices = c("All")),
                helpText("Filters both plots below to a single value of ",
                         "whichever \"Color points by\" is set to -- a ",
                         "display filter only, it doesn't change how ",
                         "TR/TC were computed."),
                hr(),
                selectInput("label_choice", "Point label:",
                            choices = c("Group" = "Group", "Rank" = "Rank"))
              ),
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
  ),

  nav_panel("Sampling Sensitivity",
            layout_sidebar(
              sidebar = sidebar(
                title = "Rarefaction settings",
                helpText("Checks whether TR and TC stay stable as fewer",
                         "taxonomic records are sampled. TR is an",
                         "unweighted mean over unique lineages, so it is",
                         "expected to stay flat. TC is abundance-weighted",
                         "and can shift if a Dataset's abundance is",
                         "concentrated in a few coarsely-identified",
                         "lineages -- this tab lets you check that",
                         "directly on your own data."),
                numericInput("rarefaction_n_boot", "Bootstrap draws per level:",
                             value = 100, min = 10, max = 1000, step = 10),
                sliderInput("rarefaction_min_frac", "Minimum sampled proportion (%):",
                            min = 5, max = 50, value = 10, step = 5),
                numericInput("rarefaction_steps", "Number of sampling levels:",
                             value = 10, min = 3, max = 20),
                actionButton("run_rarefaction", "Run rarefaction analysis",
                             class = "btn-primary w-100"),
                helpText("Datasets with fewer than 4 records are skipped.",
                         "Re-run after changing the weight vector or",
                         "editing data in the Data Editor tab -- this",
                         "analysis does not update automatically.")
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  "TR & TC vs. sampling effort",
                  downloadButton("download_rarefaction_plot", "PNG", class = "btn-sm")
                ),
                plotOutput("rarefaction_plot", height = "500px")
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  "Rarefaction data",
                  downloadButton("download_rarefaction_data", "Download CSV", class = "btn-sm")
                ),
                DTOutput("rarefaction_table")
              )
            )
  ),

  nav_panel("Rank Comparison",
            layout_sidebar(
              sidebar = sidebar(
                title = "Rank Comparison settings",
                selectInput("rank_plot_metric", "Metric (bubble size & color):",
                            choices = c("Taxonomic Resolution (TR)" = "TR",
                                        "Taxonomic Completeness (TC)" = "TC")),
                selectInput("rank_plot_rank", "Taxonomic rank (Y axis):", choices = NULL),
                numericInput("rank_plot_topn", "Show top N rank values (by abundance):",
                             value = 30, min = 3, max = 300),
                helpText("X axis follows \"Compare by\" (Data Editor ",
                         "sidebar) -- Dataset or Year. Deep ranks (Genus, ",
                         "Species) can have hundreds of distinct values; ",
                         "Top N keeps only the most abundant ones so the ",
                         "plot stays readable -- raise it to see more, or ",
                         "set it above the total count to show everything.")
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  "TR/TC by taxonomic rank",
                  downloadButton("download_rank_plot", "PNG", class = "btn-sm")
                ),
                plotOutput("rank_plot")
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  "Rank comparison data",
                  downloadButton("download_rank_data", "Download CSV", class = "btn-sm")
                ),
                DTOutput("rank_table")
              )
            )
  )
)

# --- 3. SERVER ---
server <- function(input, output, session) {
  
  # --- Example data ---
  v <- reactiveValues(data = {
    csv_text <- "Dataset;Phylum;Class;Order;Family;Genus;Species;Abundance
Dataset_Sp;P1;C1;O1;F1;G1;S1;7
Dataset_Sp;P1;C1;O2;F2;G2;S2;9
Dataset_Sp;P1;C1;O3;F3;G3;S3;500
Dataset_Sp;P1;C1;O4;F4;G4;S4;10
Dataset_Sp;P1;C1;O5;F5;G5;S5;8
Dataset_Sp;P1;C1;O6;F6;G6;S6;2
Dataset_Sp;P1;C1;O7;F7;G7;S7;100
Dataset_Class;P2;C1;;;;;100
Dataset_Abund;P3;C1;;;;;19
Dataset_Abund;P3;C1;;;;;56
Dataset_Abund;P3;C1;O3;F1;G1;S1;1200
Dataset_Abund;P3;C1;;;;;49
Dataset_Abund;P3;C1;;;;;180
Dataset_Abund;P3;C2;;;;;112
Dataset_Abund;P3;C2;;;;;85
Dataset_Abund;P3;;;;;;112
Dataset_Abund;P3;;;;;;30
Dataset_Abund;P3;;;;;;60
Dataset_Defic;P5;C1;;;;;85
Dataset_Defic;P5;C2;;;;;70
Dataset_Defic;P5;C1;O2;;;;55
Dataset_Defic;P5;C1;O2;F3;;;31
Dataset_Defic;P5;;;;;;78
Dataset_Defic;P5;;;;;;96
Dataset_Defic;P5;;;;;;91
Dataset_Defic;P5;C1;;;;;55
Dataset_Gen;P1;C1;O1;F1;G1;;70
Dataset_Gen;P1;C1;O2;F2;G2;;60
Dataset_Gen;P1;C1;O3;F3;G3;;200
Dataset_Gen;P1;C1;O4;F4;G4;;3
Dataset_Gen;P1;C1;O5;F5;G5;;70
Dataset_Ord;P5;C1;O5;;;;36
Dataset_Ord;P5;C1;O2;;;;54
Dataset_Ord;P5;C2;O1;;;;36
Dataset_Ord;P5;C1;O4;;;;120
Dataset_Ord;P5;C1;O3;;;;94
Dataset_Well_Res;P1;C1;O1;F1;G1;S1;9
Dataset_Well_Res;P1;C1;O2;F2;G2;S2;8
Dataset_Well_Res;P1;C1;O3;F3;G3;S3;3
Dataset_Well_Res;P1;C1;O4;F4;G4;S4;6
Dataset_Well_Res;P1;C1;O5;F5;G5;S5;10
Dataset_Well_Res;P1;C1;O6;F6;G6;S6;8
Dataset_Well_Res;P1;C1;O7;F7;G7;;400
Dataset_Well_Res;P1;C1;O8;F8;G8;;60
Dataset_Phy;P2;;;;;;800
Dataset_Under;P3;C1;O5;F5;;;20
Dataset_Under;P3;C1;O2;F2;G2;S2;220
Dataset_Under;P3;C1;;;;;240
Dataset_Under;P3;C1;;;;;800
Dataset_Under;P3;C1;O5;F5;G5;S5;27
Dataset_Under;P3;C1;O8;F8;G10;;10
Dataset_Under;P3;C1;O8;F8;G9;;60
Dataset_Under;P3;C1;O8;F8;G8;S8;17
Dataset_Under;P3;;;;;;200
Dataset_Under;P3;C1;O8;F8;G8;S8;1
Dataset_Fam;P5;C1;O2;F2;;;22
Dataset_Fam;P5;C1;O2;F7;;;29
Dataset_Fam;P5;C1;O3;F3;;;70
Dataset_Fam;P5;C1;O4;F4;;;69
Dataset_Fam;P5;C1;O3;F3;;;41
Dataset_Fam;P5;C1;O4;F4;;;55
Dataset_Fam;P5;C1;O5;F5;;;76
Dataset_Fam;P5;C1;O3;F3;;;64
Dataset_Fam;P5;C1;O4;F4;;;57
Dataset_Fam;P5;C1;O10;F1;;;50
Dataset_Class;P2;C2;;;;;100"
    df <- read_delim(file = csv_text, delim = ";", show_col_types = FALSE)
    as.data.frame(df)
  })
  
  # Update the initial line counting 
  observe({
    updateNumericInput(session, "row_count", value = nrow(v$data))
  })
  
  observeEvent(input$add_column, {
    new_col_name <- paste0("NewRank_", ncol(v$data) - 1)
    abund_idx <- which(names(v$data) == "Abundance")
    v$data <- v$data %>% mutate(!!new_col_name := "") %>% select(1:(abund_idx-1), last_col(), Abundance)
    updateTextInput(session, "weight_vec_str", value = paste0(input$weight_vec_str, ", ", ncol(v$data)-2))
  })

  observeEvent(input$add_year_column, {
    req(v$data)
    if (!"Year" %in% names(v$data)) {
      v$data <- v$data %>% mutate(Year = NA_character_, .after = "Dataset")
    }
  })

  current_weights <- reactive({ as.numeric(unlist(strsplit(input$weight_vec_str, ","))) })

  # "Compare by" drives every grouping/plot in the app (Summary, Biplot,
  # quadrant chart, Sampling Sensitivity, Rank Comparison's X axis).
  # "Year" only ever appears as a choice once the column exists AND has at
  # least one non-blank value -- offering a grouping option that would
  # produce nothing but an empty/NA group isn't didactic, it's confusing.
  observe({
    choices <- c("Dataset" = "Dataset")
    if ("Year" %in% names(v$data) && any(!is.na(v$data$Year) & v$data$Year != "")) {
      choices <- c(choices, "Year" = "Year")
    }
    selected <- if (!is.null(input$primary_col) && input$primary_col %in% choices) input$primary_col else "Dataset"
    updateSelectInput(session, "primary_col", choices = choices, selected = selected)
  })

  current_primary_col <- reactive({
    if (!is.null(input$primary_col) && input$primary_col %in% names(v$data)) input$primary_col else "Dataset"
  })

  observeEvent(input$upload_csv, {
    req(input$upload_csv)
    line <- readLines(input$upload_csv$datapath, n = 1)
    sep <- if(grepl(";", line)) ";" else ","
    df <- read_delim(input$upload_csv$datapath, delim = sep, show_col_types = FALSE)

    found_n <- intersect(names(df), c("individualCount", "n", "Abundance", "abundance", "count"))[1]
    if(!is.na(found_n)) df <- df %>% rename(Abundance = !!sym(found_n))
    if(!"Abundance" %in% names(df)) df$Abundance <- 1
    if(!"Dataset" %in% names(df)) df$Dataset <- "Uploaded"

    found_year <- intersect(names(df), c("Year", "year", "YEAR"))[1]
    if (!is.na(found_year) && found_year != "Year") df <- df %>% rename(Year = !!sym(found_year))

    tax_cols <- setdiff(names(df), c("Dataset", "Year", "Abundance"))
    v$data <- df %>% select(Dataset, any_of("Year"), all_of(tax_cols), Abundance) %>% as.data.frame()
    updateNumericInput(session, "row_count", value = nrow(v$data))
    updateTextInput(session, "weight_vec_str", value = paste(seq_along(tax_cols), collapse = ", "))
  })
  
  observeEvent(input$row_count, {
    req(v$data)
    curr_n <- nrow(v$data)
    if (input$row_count > curr_n) {
      extra <- v$data[rep(1, input$row_count - curr_n), ]
      extra[extra != ""] <- ""; extra$Abundance <- 1; extra$Dataset <- "New_Entry"
      v$data <- rbind(v$data, extra)
    } else if (input$row_count < curr_n) { 
      v$data <- v$data[1:input$row_count, ] 
    }
  })
  
  output$editable_table <- renderDT({
    datatable(v$data, editable = list(target = 'all'), options = list(pageLength = 66, dom = 'tip'))
  })
  
  observeEvent(input$editable_table_cell_edit, { v$data <<- editData(v$data, input$editable_table_cell_edit) })
  
  # Confidence is a display-only flag (does NOT change TR/TC): groups with
  # fewer records (Total_N, i.e. summed Abundance) than the sidebar
  # threshold are marked "Low confidence" per Reviewer #2's request about
  # sensitivity to sparse data.
  processed_results <- reactive({
    req(v$data)
    res <- calculate_TaK_shiny(v$data, current_weights(), primary_col = current_primary_col())
    req(res)
    threshold <- if (is.null(input$low_confidence_threshold)) 10 else input$low_confidence_threshold
    res$Confidence <- ifelse(res$Total_N < threshold, "Low confidence", "OK")
    res
  })

  # Taxonomic level picker for the Visualization tab specifically -- kept
  # separate from Summary/Rank Comparison's own rank choices (each tab
  # picks its own rank independently; Summary always aggregates at the
  # first taxonomic column no matter what's chosen here).
  observe({
    req(v$data)
    meta_cols <- c("Dataset", "Year", "Abundance")
    tax_cols <- setdiff(names(v$data), meta_cols)
    if (length(tax_cols) == 0) return()
    selected <- if (!is.null(input$viz_rank_col) && input$viz_rank_col %in% tax_cols)
      input$viz_rank_col else tax_cols[1]
    updateSelectInput(session, "viz_rank_col", choices = tax_cols, selected = selected)
  })

  # Same Confidence-flagging as processed_results(), but computed at
  # whatever rank_col "Taxonomic level" (Visualization tab) currently
  # picks -- feeds only the Biplot/quadrant chart, so switching it never
  # changes Summary's numbers.
  viz_results <- reactive({
    req(v$data, input$viz_rank_col)
    res <- calculate_TaK_shiny(v$data, current_weights(), primary_col = current_primary_col(),
                                rank_col = input$viz_rank_col)
    req(res)
    threshold <- if (is.null(input$low_confidence_threshold)) 10 else input$low_confidence_threshold
    res$Confidence <- ifelse(res$Total_N < threshold, "Low confidence", "OK")
    res
  })

  # Which column "Color points by" currently means -- Group (Dataset/Year)
  # or Rank (the taxonomic level above) -- shared by the Show filter, the
  # Biplot's color aesthetic and the quadrant chart's X axis, so the three
  # never disagree about what a color/bar represents.
  viz_color_col <- reactive({
    if (!is.null(input$viz_color_by) && input$viz_color_by == "Rank") "Rank" else "Group"
  })

  viz_color_label <- reactive({
    if (viz_color_col() == "Rank") input$viz_rank_col else current_primary_col()
  })

  # Local display filter for the Visualization tab -- doesn't touch how
  # TR/TC were computed, just which value(s) get plotted.
  observe({
    res <- viz_results()
    req(res)
    vals <- sort(unique(res[[viz_color_col()]]))
    selected <- if (!is.null(input$viz_group_filter) && input$viz_group_filter %in% c("All", vals))
      input$viz_group_filter else "All"
    updateSelectInput(session, "viz_group_filter", choices = c("All", vals), selected = selected)
  })

  viz_filtered_results <- reactive({
    res <- viz_results(); req(res)
    if (!is.null(input$viz_group_filter) && input$viz_group_filter != "All") {
      res <- res %>% filter(.data[[viz_color_col()]] == input$viz_group_filter)
    }
    res
  })

  summary_df <- reactive({
    res <- processed_results()
    req(res)
    out <- res %>%
      group_by(Group) %>%
      summarise(
        Mean_TR = mean(TR, na.rm = TRUE),
        Mean_TC = mean(TC, na.rm = TRUE),
        Total_Individuals = sum(Total_N, na.rm = TRUE),
        Taxa_Groups_Count = n(),
        Total_unique_taxa = sum(unique_taxa, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(across(c(Mean_TR, Mean_TC), ~round(., 3)))
    names(out)[names(out) == "Group"] <- current_primary_col()
    out
  })

  output$summary_table <- renderDT({
    datatable(summary_df(), options = list(pageLength = 10, dom = 't'))
  })

  plot_biplot <- reactive({
    res <- viz_filtered_results(); req(res)
    color_col <- viz_color_col()
    label_col <- if (!is.null(input$label_choice) && input$label_choice %in% names(res)) input$label_choice else "Group"
    ggplot(res, aes(x = TC, y = TR, color = .data[[color_col]])) +
      annotate("rect", xmin=0.5, xmax=1, ymin=0.5, ymax=1, fill="#2ecc71", alpha=0.1) +
      annotate("rect", xmin=0, xmax=0.5, ymin=0.5, ymax=1, fill="#3498db", alpha=0.1) +
      annotate("rect", xmin=0, xmax=0.5, ymin=0, ymax=0.5, fill="#e74c3c", alpha=0.1) +
      annotate("rect", xmin=0.5, xmax=1, ymin=0, ymax=0.5, fill="#f1c40f", alpha=0.1) +
      geom_point(aes(size=Total_N, shape=Confidence), alpha=0.7) +
      geom_text_repel(aes(label = .data[[label_col]])) +
      geom_abline(intercept=0, slope=1, linetype="dashed") +
      scale_shape_manual(values = c("OK" = 16, "Low confidence" = 1),
                          name = "Confidence") +
      scale_x_continuous(limits=c(0, 1.05)) + scale_y_continuous(limits=c(0, 1.05)) +
      labs(color = viz_color_label()) +
      theme_minimal() + theme(legend.position = "bottom")
  })

  plot_quadrant <- reactive({
    res <- viz_filtered_results(); req(res)
    color_col <- viz_color_col()
    quad_data <- res %>%
      mutate(Quadrant = case_when(
        TC >= 0.5 & TR >= 0.5 ~ "Well resolved",
        TC <  0.5 & TR >= 0.5 ~ "Under resolved",
        TC <  0.5 & TR <  0.5 ~ "Data deficient",
        TC >= 0.5 & TR <  0.5 ~ "Abundance biased"
      )) %>%
      group_by(.data[[color_col]], Quadrant) %>%
      summarise(Count = n(), .groups = 'drop_last') %>%
      mutate(Perc = Count / sum(Count))

    ggplot(quad_data, aes(x = .data[[color_col]], y = Perc, fill = Quadrant)) +
      geom_bar(stat = "identity", position = "fill", color = "white") +
      geom_text(aes(label = percent(Perc, accuracy = 1)), position = position_fill(vjust = 0.5), fontface="bold") +
      scale_fill_manual(values = c("Well resolved"="#2ecc71", "Under resolved"="#3498db", "Abundance biased"="#f1c40f", "Data deficient"="#e74c3c")) +
      scale_y_continuous(labels = percent) +
      labs(x = viz_color_label()) +
      theme_minimal() +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 40, hjust = 1))
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
  
  # --- Sampling Sensitivity (rarefaction) ---
  # Deliberately NOT reactive to every input change (unlike the Biplot /
  # Summary tabs) -- a bootstrap sweep over all sampling levels takes a
  # few seconds even for a mid-sized dataset, so it only reruns when the
  # user explicitly clicks "Run rarefaction analysis".
  rarefaction_fractions <- reactive({
    steps <- input$rarefaction_steps
    min_frac <- input$rarefaction_min_frac / 100
    if (is.null(steps) || is.null(min_frac) || is.na(steps) || steps < 2) {
      return(seq(0.1, 1, by = 0.1))
    }
    seq(min_frac, 1, length.out = steps)
  })

  rarefaction_result <- eventReactive(input$run_rarefaction, {
    req(v$data)
    n_boot <- if (is.null(input$rarefaction_n_boot) || is.na(input$rarefaction_n_boot))
      100 else input$rarefaction_n_boot
    fractions <- rarefaction_fractions()
    withProgress(message = "Running rarefaction analysis", value = 0, {
      rarefy_tak(v$data, current_weights(), primary_col = current_primary_col(),
                 fractions = fractions, n_boot = n_boot,
                 progress_callback = function(step, total) {
                   incProgress(1 / total, detail = sprintf("%d of %d", step, total))
                 })
    })
  })

  output$rarefaction_plot <- renderPlot({
    res <- rarefaction_result()
    validate(need(!is.null(res),
                  "No Group in the current data has 4 or more records -- rarefaction needs at least that many lineages to subsample."))
    draw_rarefaction_plot(res, group_label = current_primary_col())
  })

  output$rarefaction_table <- renderDT({
    res <- rarefaction_result()
    req(res)
    datatable(res %>% mutate(across(where(is.numeric), ~round(., 4))),
              options = list(pageLength = 10, dom = 'tip'))
  })

  output$download_rarefaction_plot <- downloadHandler(
    filename = function() { paste("rarefaction_", Sys.Date(), ".png", sep = "") },
    content = function(file) {
      res <- rarefaction_result(); req(res)
      ggsave(file, plot = draw_rarefaction_plot(res, group_label = current_primary_col()),
             width = 10, height = 5.5, dpi = 300)
    }
  )

  output$download_rarefaction_data <- downloadHandler(
    filename = function() { paste("rarefaction_data_", Sys.Date(), ".csv", sep = "") },
    content = function(file) {
      res <- rarefaction_result(); req(res)
      write.csv(res, file, row.names = FALSE)
    }
  )

  # --- Rank Comparison ---
  observe({
    req(v$data)
    meta_cols <- c("Dataset", "Year", "Abundance")
    tax_cols <- setdiff(names(v$data), meta_cols)
    if (length(tax_cols) == 0) return()
    selected <- if (!is.null(input$rank_plot_rank) && input$rank_plot_rank %in% tax_cols)
      input$rank_plot_rank else tax_cols[1]
    updateSelectInput(session, "rank_plot_rank", choices = tax_cols, selected = selected)
  })

  rank_plot_result <- reactive({
    req(v$data, input$rank_plot_rank)
    calculate_TaK_shiny(v$data, current_weights(), primary_col = current_primary_col(),
                         rank_col = input$rank_plot_rank)
  })

  rank_plot_topn <- reactive({
    n <- input$rank_plot_topn
    if (is.null(n) || is.na(n) || n <= 0) 30 else n
  })

  plot_rank_comparison <- reactive({
    res <- rank_plot_result(); req(res)
    metric <- if (is.null(input$rank_plot_metric)) "TR" else input$rank_plot_metric
    draw_rank_comparison_plot(res, metric = metric, group_label = current_primary_col(),
                               rank_label = input$rank_plot_rank, top_n = rank_plot_topn())
  })

  # A fixed plot height made rows bleed into each other once a rank had
  # more than ~15-20 distinct values (Genus/Species easily have hundreds) --
  # height now scales with how many rows are actually shown (after the Top
  # N cap), on screen and in the PNG download alike, via the same
  # rank_comparison_height_*() helpers so the two never drift apart.
  rank_plot_n_shown <- reactive({
    res <- rank_plot_result(); req(res)
    rank_comparison_n_shown(res, rank_plot_topn())
  })

  output$rank_plot <- renderPlot({
    plot_rank_comparison()
  }, height = function() rank_comparison_height_px(rank_plot_n_shown()))

  output$rank_table <- renderDT({
    res <- rank_plot_result(); req(res)
    datatable(res %>% mutate(across(c(TR, TC), ~round(., 4))),
              options = list(pageLength = 10, dom = 'tip'))
  })

  output$download_rank_plot <- downloadHandler(
    filename = function() { paste("rank_comparison_", Sys.Date(), ".png", sep = "") },
    content = function(file) {
      ggsave(file, plot = plot_rank_comparison(), width = 10,
             height = rank_comparison_height_in(rank_plot_n_shown()), dpi = 300, bg = "white")
    }
  )

  output$download_rank_data <- downloadHandler(
    filename = function() { paste("rank_comparison_data_", Sys.Date(), ".csv", sep = "") },
    content = function(file) {
      res <- rank_plot_result(); req(res)
      write.csv(res, file, row.names = FALSE)
    }
  )

  observeEvent(input$reset_data, { session$reload() })
}

shinyApp(ui, server)
