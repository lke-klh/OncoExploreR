library(shiny)
library(tidyverse)
library(ggplot2)
library(plotly)
library(pheatmap)
library(survival)
library(randomForestSRC)

simulate_gene_effect <- function(
    rsf_fit,
    clean_df, 
    age_input,
    gender_input,
    expr_baseline,
    expr_treated,
    N = 500,
    years = c(3, 5, 10),
    expr_sd = 5.0
) {
  
  gender_factor <- factor(gender_input, levels = levels(clean_df$gender))
  time_grid <- rsf_fit$time.interest
  gene_range <- range(clean_df$Gene_Expression, na.rm = TRUE)
  
  clamp_gene <- function(x) {
    pmin(pmax(x, gene_range[1]), gene_range[2])
  }
  
  new_patients <- data.frame(
    age = c(age_input, age_input),
    gender = gender_factor,
    Gene_Expression = c(expr_baseline, expr_treated)
  )
  
  pred_single <- predict(rsf_fit, newdata = new_patients)
  surv_mat_single <- pred_single$survival
  
  single_curve_df <- data.frame(
    time = rep(time_grid, times = nrow(surv_mat_single)),
    survival = as.vector(t(surv_mat_single)),
    group = factor(
      rep(c("Baseline expression", "Treated expression"),
          each = length(time_grid)),
      levels = c("Baseline expression", "Treated expression")
    )
  )
  
  # Simulation
  gene_base_vec <- clamp_gene(rnorm(N, mean = expr_baseline, sd = expr_sd))
  gene_treat_vec <- clamp_gene(rnorm(N, mean = expr_treated,  sd = expr_sd))
  
  base_patients <- data.frame(
    age = rep(age_input, N),
    gender = rep(gender_factor, N),
    Gene_Expression = gene_base_vec
  )
  
  treat_patients <- data.frame(
    age = rep(age_input, N),
    gender = rep(gender_factor, N),
    Gene_Expression = gene_treat_vec
  )
  
  pred_base <- predict(rsf_fit, newdata = base_patients)
  pred_treat <- predict(rsf_fit, newdata = treat_patients)
  
  S_base <- pred_base$survival
  S_treat <- pred_treat$survival
  
  mean_base <- colMeans(S_base)
  mean_treat <- colMeans(S_treat)
  
  q25_base <- apply(S_base,  2, quantile, 0.25)
  q75_base <- apply(S_base,  2, quantile, 0.75)
  q25_treat <- apply(S_treat, 2, quantile, 0.25)
  q75_treat <- apply(S_treat, 2, quantile, 0.75)
  
  avg_curve_df <- rbind(
    data.frame(
      time = time_grid,
      mean = mean_base,
      lower = q25_base,
      upper = q75_base,
      scenario = "Baseline"
    ),
    data.frame(
      time = time_grid,
      mean = mean_treat,
      lower = q25_treat,
      upper = q75_treat,
      scenario = "Treated"
    )
  )
  
  # survival gain
  t_star <- years * 365
  idx_vec <- sapply(t_star, function(t) which.min(abs(time_grid - t)))
  
  S_base_years <- S_base[, idx_vec, drop = FALSE]
  S_treat_years <- S_treat[, idx_vec, drop = FALSE]
  
  delta_S <- S_treat_years - S_base_years
  
  delta_hist_df <- data.frame(
    delta = as.vector(delta_S),
    year = factor(rep(years, each = N))
  )
  
  summary_list <- lapply(seq_along(years), function(j) {
    list(
      year = years[j],
      mean_diff = mean(delta_S[, j]),
      quantiles = quantile(delta_S[, j], c(0.1, 0.5, 0.9))
    )
  })
  
  list(
    single_curve_df = single_curve_df,
    avg_curve_df = avg_curve_df,
    delta_hist_df = delta_hist_df,
    summary = summary_list
  )
}


server <- function(input, output, session) {

  classic_purple_gray <- c(
    "#7C75C4", "#C9A4BA", "#E3B4C2",
    "#EBD4D0", "#7E7C48", "#A39F53",
    "#BFA7D5", "#D0B7E1", "#CDBDBC"
  )
  
  output$organ_title <- renderText({
    req(input$selected_organ)
    paste(input$selected_organ)
  })
  
  output$organ_intro <- renderUI({
    req(input$selected_organ)
    text <- switch(
      input$selected_organ,
      "Thyroid" = tagList(
        p(strong(em("Thyroid cancer"))),
        p("Thyroid cancer arises from thyroid cells in the neck. 
          Most tumors are slow-growing and highly curable, often found as a painless nodule. 
          Prognosis is generally excellent, especially for younger patients and early-stage disease.")
      ),
      
      "Bronchus and Lung" = tagList(
        p(strong(em("Lung cancer"))),
        p("Lung cancer develops in the bronchi and lung tissue. 
          Smoking is the main risk factor, but air pollution and genetics also contribute. 
          It often presents with persistent cough, shortness of breath, chest pain, or weight loss.")
      ),
      
      "Breast" = tagList(
        p(strong(em("Breast cancer"))),
        p("Breast cancer begins in breast ducts or lobules. 
          It is common, especially among women, and may appear as a lump, skin change, or nipple discharge. 
          Hormone receptor and HER2 status guide targeted and systemic therapies.")
      ),
      
      "Liver" = tagList(
        p(strong(em("Primary liver cancer"))),
        p("Primary liver cancer, often hepatocellular carcinoma, 
        usually arises in chronically damaged livers from hepatitis B, hepatitis C, or cirrhosis. 
          Symptoms can be subtle, such as fatigue, abdominal discomfort, or weight loss, 
          making surveillance in high-risk groups important.")
      ),
      
      "Kidney" = tagList(
        p(strong(em("Kidney cancer"))),
        p("Kidney cancer, typically renal cell carcinoma, forms in the renal cortex. 
          It may be discovered incidentally or present with blood in urine, flank pain, or a mass. 
          Smoking, obesity, hypertension, and some hereditary syndromes increase risk.")
      ),
      
      "Colon" = tagList(
        p(strong(em("Colorectal cancer"))),
        p("Colorectal cancer usually develops from polyps in the colon or rectum over years. 
          Screening colonoscopy can detect and remove precancerous lesions. 
          Symptoms include blood in stool, altered bowel habits, anemia, or abdominal pain, 
          though early disease may be asymptomatic.")
      )
    )
  })
  
  # Gene Heat Map
  output$plot_heatmap <- renderPlot({
    req(input$selected_organ)
    expr_data <- get_cancer_expr(input$selected_organ)
      
    n_genes <- min(20, nrow(expr_data))
    top_genes_data <- expr_data[1:n_genes, , drop = FALSE]
    
    if (ncol(top_genes_data) > 50) {
      set.seed(123)
      keep_idx <- sample(seq_len(ncol(top_genes_data)), 50)
      top_genes_data <- top_genes_data[, keep_idx, drop = FALSE]
    }
      
    top_genes_data[is.na(top_genes_data)] <- 0
    top_genes_data[is.infinite(top_genes_data)] <- 0
      
    pheatmap(
      top_genes_data,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      scale = "none",
      color = colorRampPalette(c("steelblue", "white", "darkred"))(50),
      main = paste0("Top ", n_genes, " Variable Genes (from 50 Random Tumor Samples)"),
      show_colnames = FALSE,
      show_rownames = TRUE,
      fontsize_row = 10
    )
  })
  
  
  # Race
  output$plot_race <- renderPlotly({
    req(input$selected_organ)
    df <- get_cancer_data(input$selected_organ)
    req(df)
    
    race_df <- df %>%
      count(race) %>%
      mutate(race = reorder(race, n, decreasing = TRUE))
    
    plot_ly(race_df, 
            x = ~race, 
            y = ~n, 
            type = "bar",
            color = ~race, 
            colors = classic_purple_gray, 
            text = ~n,
            textposition = "outside",
            hoverinfo = "x+y") %>%
      layout(
        title = list(text = "Race Distribution", font = list(size = 14)),
        xaxis = list(title = "", tickangle = -15),
        yaxis = list(title = "Count"),
        margin = list(t = 60, b = 60, l = 60, r = 20),
        showlegend = FALSE
      )
  })
  
  # Sex 
  output$plot_gender <- renderPlotly({
    req(input$selected_organ)
    df <- get_cancer_data(input$selected_organ)
    req(df)
    
    gender_df <- df %>%
      count(gender) %>%
      mutate(
        perc = n / sum(n),
        label = paste0(round(perc * 100, 1), "%")
      )
    
    plot_ly(gender_df, 
            labels = ~gender, 
            values = ~n, 
            type = "pie",
            textposition = "inside",
            textinfo = "label+percent",
            marker = list(colors = classic_purple_gray),
            sort = FALSE) %>%
      layout(
        title = list(text = "Sex Distribution", font = list(size = 14)),
        margin = list(t = 50, b = 20, l = 20, r = 20),
        showlegend = TRUE
      )
  })
  
  # Age
  output$plot_age <- renderPlotly({
    req(input$selected_organ)
    df <- get_cancer_data(input$selected_organ)
    req(df)
    
    breaks <- c(0, 50, 60, 70, 80, 120)
    labels <- c("<50", "50-59", "60-69", "70-79", "80+")
    
    age_df <- df %>%
      filter(!is.na(age)) %>%
      mutate(
        age_group = cut(
          age,
          breaks = breaks,
          labels = labels,
          right = FALSE
        )
      ) %>%
      count(age_group)
    
    plot_ly(age_df,
            x = ~age_group,
            y = ~n,
            type = "bar",
            color = ~age_group,
            colors = classic_purple_gray,
            text = ~n,
            textposition = "outside",
            hoverinfo = "x+y") %>%
      layout(
        title = list(text = "Age Distribution", font = list(size = 14)),
        xaxis = list(title = "Age Group"),
        yaxis = list(title = "Count"),
        margin = list(t = 60, b = 60, l = 60, r = 20),
        showlegend = FALSE
      )
  })
  
  # Survival Analysis
  sim_ready <- reactiveVal(FALSE)
  
  clean_df_surv_event <- eventReactive(input$run_sim, {
    cancer_type <- isolate(input$cancer_type_surv)
    gene <- isolate(input$selected_gene)
    
    obj <- get_cancer_object(cancer_type)
    req(obj)
    
    df <- obj$merged_data
    req(df)
    
    req(all(c("tissue_status", "patient_id", "time",
              "event", "age", "gender") %in% colnames(df)))
    req(gene %in% colnames(df))
    
    df %>%
      filter(tissue_status == "Primary Tumor") %>%
      dplyr::select(
        patient_id, time, event, age, gender,
        Gene_Expression = all_of(gene)
      ) %>%
      mutate(
        event = dplyr::case_when(
          event == "Dead"  ~ 1,
          event == "Alive" ~ 0,
          TRUE ~ NA_real_
        ),
        gender = as.factor(gender),
        time = as.numeric(time),
        age = as.numeric(age)
      ) %>%
      na.omit()
  })
  
  rsf_fit_event <- eventReactive(input$run_sim, {
    clean_df <- clean_df_surv_event()
    set.seed(123)
    rfsrc(
      Surv(time, event) ~ age + gender + Gene_Expression,
      data = clean_df,
      ntree = 200,
      importance = TRUE
    )
  })
  
  sim_res_event <- eventReactive(input$run_sim, {
    rsf_fit <- rsf_fit_event()
    clean_df <- clean_df_surv_event()
    
    gender_input <- isolate(input$gender_input)
    age_input <- isolate(input$age_input)
    expr_treated <- isolate(input$expr_treated)
    
    expr_baseline <- median(
      clean_df$Gene_Expression[clean_df$gender == gender_input],
      na.rm = TRUE
    )
    if (is.na(expr_baseline)) {
      expr_baseline <- median(clean_df$Gene_Expression, na.rm = TRUE)
    }
    
    res <- simulate_gene_effect(
      rsf_fit = rsf_fit,
      clean_df = clean_df,
      age_input = age_input,
      gender_input = gender_input,
      expr_baseline = expr_baseline,
      expr_treated = expr_treated,
      N = 500,
      years = c(3, 5, 10),
      expr_sd = 5.0
    )
    
    sim_ready(TRUE)
    res
  })
  

  output$surv_curv_ui <- renderUI({
    if (input$run_sim == 0) {
      div(
        class = "panel-body",
        p("Select a cancer type, sex, age, gene, and expression level, then click "),
        p(strong("Run Survival Analysis")),
        p("to generate survival curves.")
      )
    } else {
      tagList(
        h3(
          span("Survival Curves for Gene ", style = "font-weight: 400;"),
          span(input$selected_gene, style = "color: #911D2A; font-weight: 400;"),
          span(" in ", style = "font-weight: 400;"),
          span(input$cancer_type_surv, style = "color: #911D2A; font-weight: 400;"),
          span("Cancer", style = "font-weight: 400;")
        ),
        withSpinner(plotlyOutput("surv_curv", height = "500px"),
                    type = 4,
                    color = "#7A658A",
                    size = 1)
      )
    }
  })
  
  output$surv_gain_ui <- renderUI({
    if (input$run_sim == 0) {
      div(
        class = "panel-body",
        p("Run the simulation to see survival gains at 3, 5, and 10 years.")
      )
    } else {
      tagList(
        h3(
          span("Survival Gains for Gene ", style = "font-weight: 400;"),
          span(input$selected_gene, style = "color: #911D2A; font-weight: 400;"),
          span(" in ", style = "font-weight: 400;"),
          span(input$cancer_type_surv, style = "color: #911D2A; font-weight: 400;"),
          span("Cancer", style = "font-weight: 400;")
        ),
        withSpinner(plotlyOutput("surv_gain", height = "500px"),
                    type = 4,
                    color = "#7A658A",
                    size = 1)
      )
    }
  })
  
  
  # Plots
  output$surv_curv <- renderPlotly({
    res <- sim_res_event()
    req(res)
    
    scenario_cols <- c(
      "Baseline" = "#CC6869",
      "Treated" = "#68B8CC"
    )
    
    p <- ggplot(res$avg_curve_df,
      aes(x = time, y = mean, colour = scenario, fill = scenario)) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.35, show.legend = FALSE) +
      geom_line(linewidth = 1.2) +
      scale_colour_manual(values = scenario_cols) +
      scale_fill_manual(values = scenario_cols) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Time (days)", y = "Survival probability") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$surv_gain <- renderPlotly({
    res <- sim_res_event()
    req(res)
    
    plot_df <- res$delta_hist_df %>%
      mutate(year_label = paste0(year, " Years"),
             year_label = fct_reorder(year_label, as.numeric(year)))
    
    p <- ggplot(plot_df, aes(x = delta)) +
      geom_histogram(bins = 30, alpha = 0.7, fill = "#696175") +
      geom_vline(xintercept = 0,
                 linetype = "dashed", colour = "#CF303B") +
      facet_wrap(~ year_label, nrow = 1) +
      labs(
        x = "Change in Survival Probability (Treated vs. Baseline)",
        y = "Count"
      ) +
      theme_minimal() +
      theme(plot.margin = margin(t = 10, r = 10, b = 20, l = 10),
            strip.text = element_text(size = 12, face = "bold"))
    
    ggplotly(p)
  })

  
  # Download Data
  merged_data_reactive <- reactive({
    req(input$cancer_type_dl) 
    obj <- get_cancer_object(input$cancer_type_dl)
    return(obj$merged_data)
  })
  
  output$dl_data_preview <- renderTable({
    df <- merged_data_reactive()
    head(df, 20)
  })
  
  output$download_data_btn <- downloadHandler(
    filename = function() {
      paste0(input$cancer_type_dl, "_Merged_Data.csv")
    },
    content = function(file) {
      write.csv(merged_data_reactive(), file, row.names = FALSE)
    }
  )
}
