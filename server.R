library(shiny)
library(tidyverse)
library(ggplot2)
library(plotly)
library(pheatmap)
library(survival)
library(randomForestSRC)
library(edgeR)
library(limma)
library(DESeq2)
library(pROC)

deg_simulation <- function(data, N, pct,
                           logFC_vec = c(0.5, 1.0, 1.5),
                           phi = 0.2,
                           seed = 1234) {
  set.seed(seed)
  
  counts <- round(pmax(2^as.matrix(data) - 1, 0))
  G <- nrow(counts)
  if (G == 0) stop("No genes in input data.")
  
  genes <- rownames(counts)
  base_mu <- pmax(rowMeans(counts), 0.1)
  size_g <- 1 / phi
  
  # Differentially expressed genes
  n_DE <- round(G * pct)
  DE_idx <- sample(G, n_DE)
  DE_up  <- DE_idx[seq_len(ceiling(n_DE / 2))]
  DE_down <- setdiff(DE_idx, DE_up)
  
  direction <- rep(0, G)
  direction[DE_up]   <- 1
  direction[DE_down] <- -1
  names(direction) <- genes
  
  group <- factor(rep(c("healthy", "patient"), each = N))
  sample_names <- c(paste0("healthy", seq_len(N)),
                    paste0("patient", seq_len(N)))
  
  simulate_condition <- function(mu) {
    matrix(
      rnbinom(G * N, mu = rep(mu, each = N), size = size_g),
      nrow = G, ncol = N, byrow = TRUE
    )
  }
  
  expr_data <- lapply(logFC_vec, function(lfc) {
    fc <- 2^lfc
    
    muA <- base_mu
    muB <- base_mu
    muB[DE_up]   <- muA[DE_up] * fc
    muB[DE_down] <- muA[DE_down] / fc
    
    sim_A <- simulate_condition(muA)
    sim_B <- simulate_condition(muB)
    
    sim <- cbind(sim_A, sim_B)
    rownames(sim) <- genes
    colnames(sim) <- sample_names
    
    list(counts = sim, group = group, samples = sample_names)
  })
  
  names(expr_data) <- paste0("logFC_", logFC_vec)
  
  list(
    expr_data = expr_data,
    DE_list = genes[DE_idx],
    DE_up = genes[DE_up],
    DE_down = genes[DE_down],
    direction = direction
  )
}

simulate_gene_effect <- function(
    rsf_fit, clean_df,
    age_input, gender_input,
    expr_baseline, expr_treated,
    N = 500,
    years = c(3, 5, 10),
    expr_sd = 5
) {
  gender_factor <- factor(gender_input, levels = levels(clean_df$gender))
  time_grid <- rsf_fit$time.interest
  if (length(time_grid) == 0) stop("RSF has no time grid")
  
  # Clamp gene values to observed range
  clamp <- function(x) pmin(pmax(x, min(clean_df$Gene_Expression)), 
                            max(clean_df$Gene_Expression))

# Predict survival
surv_pred <- function(df) predict(rsf_fit, newdata = df)$survival

df_single <- data.frame(
  age = rep(age_input, 2),
  gender = gender_factor,
  Gene_Expression = c(expr_baseline, expr_treated)
)

S_single <- surv_pred(df_single)

single_curve_df <- data.frame(
  time = rep(time_grid, each = 2),
  survival = as.vector(t(S_single)),
  group = factor(rep(c("Baseline expression", "Treated expression"),
                     times = length(time_grid)))
)

# simulation
sim_patients <- function(expr_mean) {
  data.frame(
    age = rep(age_input, N),
    gender = gender_factor,
    Gene_Expression = clamp(rnorm(N, mean = expr_mean, sd = expr_sd))
  )
}

S_base <- surv_pred(sim_patients(expr_baseline))
S_treat <- surv_pred(sim_patients(expr_treated))

# survival curves
summarize_survival <- function(S, scenario) {
  data.frame(
    time = time_grid,
    mean  = colMeans(S),
    lower = apply(S, 2, quantile, 0.25),
    upper = apply(S, 2, quantile, 0.75),
    scenario = scenario
  )
}

avg_curve_df <- rbind(
  summarize_survival(S_base, "Baseline"),
  summarize_survival(S_treat, "Treated")
)

# survival gains
idx_years <- sapply(years * 365, function(t) which.min(abs(time_grid - t)))
delta <- S_treat[, idx_years, drop = FALSE] - 
  S_base[, idx_years, drop = FALSE]

delta_hist_df <- data.frame(
  delta = as.vector(delta),
  year = factor(rep(years, each = N))
)

summary_list <- lapply(seq_along(years), function(i) {
  list(
    year = years[i],
    mean_diff = mean(delta[, i]),
    quantiles = quantile(delta[, i], c(0.1, 0.5, 0.9))
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
  
  bodymap_data <- reactive({
    req(input$selected_organ)
    get_merged_data(input$selected_organ)
  })
  
  organ_ready <- eventReactive(input$selected_organ, {
    df <- bodymap_data()
    df$race <- as.character(df$race)
    df$gender <- as.character(df$gender)
    df
  })
  
  
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
    expr_data <- get_expr_log(input$selected_organ)
    
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
      color = colorRampPalette(c("#2553A8", "white", "#A8253F"))(50),
      main = paste0("Top ", n_genes, " Variable Genes (from 50 Random Tumor Samples)"),
      show_colnames = FALSE,
      show_rownames = TRUE,
      fontsize_row = 10
    )
  })
  
  
  # Race
  output$plot_race <- renderPlotly({
    req(input$analysis_tabs == "Race Breakdown")
    df <- organ_ready()
    req(df)
    
    race_df <- df %>%
      dplyr::group_by(race) %>%
      dplyr::summarize(n = n()) %>% 
      dplyr::mutate(race = reorder(race, n, decreasing = TRUE))
    
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
    req(input$analysis_tabs == "Sex Breakdown")
    df <- organ_ready()
    req(df)
    
    gender_df <- df %>%
      dplyr::count(gender) %>%
      dplyr::mutate(
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
    req(input$analysis_tabs == "Age Breakdown")
    df <- organ_ready()
    req(df)
    
    breaks <- c(0, 50, 60, 70, 80, 120)
    labels <- c("<50", "50-59", "60-69", "70-79", "80+")
    
    age_df <- df %>%
      dplyr::filter(!is.na(age)) %>%
      dplyr::mutate(
        age_group = cut(
          age,
          breaks = breaks,
          labels = labels,
          right = FALSE
        )
      ) %>%
      dplyr::count(age_group)
    
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
  
  # Benchmark
  output$de_plots_ui <- renderUI({
    if (input$run_de == 0) {
      div(
        class = "panel-body",
        p("Select a tumor site from the sidebar, 
          adjust the sample size, signal strength, and dispersion, then click"),
        p(strong("Run DE Genes Detection")),
        p("to start the simulation.")
      )
    } else {
      tagList(
        h3("Benchmark Results"),
        tabsetPanel(
          id = "benchmark_tabs",
          
          # Tab 1: ROC Curve
          tabPanel(
            title = "ROC Curve",
            br(),
            withSpinner(plotlyOutput("bench_plot_roc", height = "550px"), 
                        type = 4, color = "#7A658A", size = 1)
          ),
          
          # Tab 2: FDR vs Sensitivity
          tabPanel(
            title = "FDR vs Sensitivity",
            br(),
            withSpinner(plotlyOutput("bench_plot_fdr", height = "550px"), 
                        type = 4, color = "#7A658A", size = 1)
          ),
          
          # Tab 3: Sensitivity vs LFC
          tabPanel(
            title = "Sensitivity vs LFC",
            br(),
            withSpinner(plotlyOutput("bench_plot_lfc", height = "550px"), 
                        type = 4, color = "#7A658A", size = 1)
          )
        )
      )
    }
  })
  
  method_cols <- c(
    "edgeR" = "#7E7C48", 
    "DESeq2" = "#BF685A", 
    "limma-voom" = "#7C75C4", 
    "limma-trend" = "#E3B4C2"
  )
  
  # ROC Curve (sensitivity vs. 1-specificity)
  output$bench_plot_roc <- renderPlotly({
    res_obj <- benchmark_res()
    req(res_obj)
    
    df <- res_obj$result_df
    
    # Calculate ROC per method
    roc_df <- df %>%
      group_by(method) %>%
      do({
        roc_obj <- roc(.$truth, -.$pvalue, quiet = TRUE)
        data.frame(
          method = unique(.$method),
          tpr = roc_obj$sensitivities,
          fpr = 1 - roc_obj$specificities
        )
      }) %>%
      ungroup()
    
    p <- ggplot(roc_df, aes(fpr, tpr, color = method)) +
      geom_line(linewidth = 1) +
      geom_abline(linetype = "dashed", color = "#9990A6") +
      scale_color_manual(values = method_cols, drop = FALSE) +
      labs(title = "ROC Curve", 
           x = "1 - Specificity (FPR)", 
           y = "Sensitivity (TPR)") +
      theme_bw()
    
    ggplotly(p)
  })
  
  # Sensitivity vs FDR
  output$bench_plot_fdr <- renderPlotly({
    res_obj <- benchmark_res()
    req(res_obj)
    
    df <- res_obj$result_df
    alpha_vec <- seq(0.01, 0.1, length.out = 10)
    
    # Calculate stats across different alpha thresholds
    summ <- df %>%
      tidyr::expand_grid(alpha = alpha_vec) %>%
      group_by(method, LFC, alpha) %>%
      summarise(
        TP = sum(truth & (padj < alpha), na.rm = TRUE),
        FP = sum(!truth & (padj < alpha), na.rm = TRUE),
        n_true = sum(truth, na.rm = TRUE),
        n_call = sum(padj < alpha, na.rm = TRUE),
        sens = ifelse(n_true == 0, NA_real_, TP / n_true),
        fdr  = ifelse(n_call == 0, NA_real_, FP / (TP + FP)),
        .groups = "drop"
      )
    
    # filter for LFC = 1
    summ_lfc1 <- summ %>% 
      mutate(LFC = suppressWarnings(as.numeric(LFC))) %>% 
      filter(LFC == 1)
    
    p <- ggplot(summ_lfc1, aes(fdr, sens, color = method, shape = method)) +
      geom_point(size = 3) +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "#9990A6") +
      scale_color_manual(values = method_cols, drop = FALSE) +
      labs(title = "Sensitivity vs FDR (LFC = 1.0)",
           x = "Observed FDR", y = "Sensitivity") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # Sensitivity vs LFC
  output$bench_plot_lfc <- renderPlotly({
    res_obj <- benchmark_res()
    req(res_obj)
    
    df <- res_obj$result_df
    
    summ_lfc <- df %>%
      group_by(method, LFC) %>%
      summarise(
        TP = sum(truth & (padj < 0.05), na.rm = TRUE),
        n_true = sum(truth, na.rm = TRUE),
        sens = ifelse(n_true == 0, 0, TP / n_true),
        .groups = "drop"
      )
    
    p <- ggplot(summ_lfc, aes(x = LFC, y = sens, color = method, marker = method)) +
      geom_point(size = 3) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = method_cols, drop = FALSE) +
      scale_x_continuous(breaks = sort(unique(summ_lfc$LFC))) +
      labs(
        title = "Sensitivity vs LFC (FDR < 0.05)",
        x = "Log Fold Change (LFC)", 
        y = "Sensitivity"
      ) +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # Sensitivity across different fold changes
  benchmark_res <- eventReactive(input$run_de, {
    req(input$cancer_type_de)
    
    N_sim <- input$de_sample_size 
    pct_sim <- input$de_pct        
    phi_sim <- input$dispersion
    
    bg_data <- get_expr_log(input$cancer_type_de)
    bg_data <- bg_data[1:500, , drop = FALSE]
    
    withProgress(message = "Running Benchmark...", value = 0, {
      incProgress(0.1, detail = "Simulating Data")
      
      # Run Simulation
      sim_obj <- deg_simulation(
        data = bg_data,
        N = N_sim,
        pct = pct_sim / 100, 
        phi = phi_sim
      )
      
      incProgress(0.4, detail = "Running DE Methods")
      
      out_df_list <- list()
      
      for (scen_name in names(sim_obj$expr_data)) {
        
        curr_sim <- sim_obj$expr_data[[scen_name]]
        counts <- curr_sim$counts
        group  <- curr_sim$group
        current_lfc <- as.numeric(gsub("logFC_", "", scen_name))
        
        design <- model.matrix(~group)
        
        # edgeR
        dge <- DGEList(counts = counts, group = group)
        dge <- calcNormFactors(dge)
        dge <- estimateDisp(dge, design)
        fit <- glmQLFit(dge, design)
        qlf <- glmQLFTest(fit, coef = 2)
        res_edger <- topTags(qlf, n = nrow(counts))$table
        
        df_edger <- data.frame(
          gene = rownames(res_edger),
          pvalue = res_edger$PValue,
          padj = res_edger$FDR,
          method = "edgeR",
          LFC = current_lfc
        )
        
        # DESeq2
        dds <- DESeqDataSetFromMatrix(countData = counts, 
                                      colData = data.frame(group=group), 
                                      design = ~group)
        
        dds <- estimateSizeFactors(dds, type = "poscounts")
        dds <- estimateDispersions(dds, quiet = TRUE)
        dds <- nbinomWaldTest(dds, quiet = TRUE)
        res_deseq <- results(dds)
        
        df_deseq <- data.frame(
          gene = rownames(res_deseq),
          pvalue = res_deseq$pvalue,
          padj = res_deseq$padj,
          method = "DESeq2",
          LFC = current_lfc
        )
        dge_lim <- DGEList(counts = counts, group = group)
        dge_lim <- calcNormFactors(dge_lim)
        
        # limma-voom
        v <- voom(dge_lim, design, plot = FALSE)
        fit_voom <- lmFit(v, design)
        fit_voom <- eBayes(fit_voom)
        res_voom <- topTable(fit_voom, coef = 2, number = nrow(counts), sort.by = "none")
        
        df_voom <- data.frame(
          gene = rownames(res_voom),
          pvalue = res_voom$P.Value,
          padj = res_voom$adj.P.Val,
          method = "limma-voom",
          LFC = current_lfc
        )
        
        # limma-trend 
        logCPM <- cpm(dge_lim, log=TRUE, prior.count=3)
        fit_trend <- lmFit(logCPM, design)
        fit_trend <- eBayes(fit_trend, trend=TRUE)
        res_trend <- topTable(fit_trend, coef = 2, number = nrow(counts), sort.by = "none")
        
        df_trend <- data.frame(
          gene = rownames(res_trend),
          pvalue = res_trend$P.Value,
          padj = res_trend$adj.P.Val,
          method = "limma-trend",
          LFC = current_lfc
        )
        out_df_list[[scen_name]] <- rbind(df_edger, df_deseq, df_voom, df_trend)
      }
      
      result_df <- do.call(rbind, out_df_list)
      
      result_df$truth <- result_df$gene %in% sim_obj$DE_list
      result_df$pvalue[is.na(result_df$pvalue)] <- 1
      result_df$padj[is.na(result_df$padj)] <- 1
      
      incProgress(0.9, detail = "Analysis Complete")
    })
    
    list(
      result_df = result_df, 
      sim_details = sim_obj
    )
  }) 
  
  # Survival Analysis
  sim_ready <- reactiveVal(FALSE)
  
  observeEvent(input$cancer_type_surv, {
    valid_genes <- get_gene_list(input$cancer_type_surv)
    
    updateSelectizeInput(
      session = session,
      inputId = "selected_gene",
      choices = valid_genes,
      selected = valid_genes[1],
      server = TRUE
    )
  })
  
  surv_analysis <- eventReactive(input$run_sim, {
    withProgress(message = "Running Survival Analysis...", value = 0, {
      
      incProgress(0.2, detail = "Loading and cleaning data")
      clean_df <- get_surv_df(input$cancer_type_surv)
      gene_col <- input$selected_gene
      clean_df$Gene_Expression <- suppressWarnings(as.numeric(clean_df[[gene_col]]))
      
      incProgress(0.45, detail = "Fitting Random Survival Forest")
      # fit RSF
      rsf_fit <- rfsrc(
        Surv(time, event) ~ age + gender + Gene_Expression,
        data = clean_df,
        ntree = 200,
        nodesize = 15,
        importance = FALSE
      )
      
      incProgress(0.75, detail = "Simulating expression scenarios")
      # run simulation
      out <- simulate_gene_effect(
        rsf_fit = rsf_fit,
        clean_df = clean_df,
        age_input = input$age_input,
        gender_input = input$gender_input,
        expr_baseline = median(clean_df$Gene_Expression),
        expr_treated = input$expr_treated,
        N = 200,               
        years = c(3, 5, 10),
        expr_sd = 5
      )
      
      incProgress(1, detail = "Finalizing results")
      # free memory
      rm(rsf_fit, clean_df); gc()
      
      return(out)
    })
  })
  
  surv_title_gene <- eventReactive(input$run_sim, {
    isolate(input$selected_gene)
  })
  
  surv_title_cancer <- eventReactive(input$run_sim, {
    isolate(input$cancer_type_surv)
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
          span(surv_title_gene(), style = "color: #911D2A; font-weight: 400;"),
          span(" in ", style = "font-weight: 400;"),
          span(surv_title_cancer(), style = "color: #911D2A; font-weight: 400;"),
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
          span(surv_title_gene(), style = "color: #911D2A; font-weight: 400;"),
          span(" in ", style = "font-weight: 400;"),
          span(surv_title_cancer(), style = "color: #911D2A; font-weight: 400;"),
          span("Cancer", style = "font-weight: 400;")
        ),
        withSpinner(plotlyOutput("surv_gain", height = "400px"),
                    type = 4,
                    color = "#7A658A",
                    size = 1)
      )
    }
  })
  
  
  # Plots
  output$surv_curv <- renderPlotly({
    res <- surv_analysis()
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
    res <- surv_analysis()
    req(res)
    
    plot_df <- res$delta_hist_df %>%
      mutate(year_label = paste0(year, " Years"),
             year_label = fct_reorder(year_label, as.numeric(year)))
    
    p <- ggplot(plot_df, aes(x = delta)) +
      geom_histogram(bins = 30, alpha = 0.7, fill = "#9990A6") +
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
    get_merged_data(input$cancer_type_dl)
  })
  
  output$dl_data_preview <- renderTable({
    withProgress(message = "Loading preview...", value = 1, {
      df <- merged_data_reactive()
      df[1:20, 1:31]
    })
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