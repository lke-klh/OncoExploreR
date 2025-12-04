library(shiny)
library(tidyverse)
library(ggplot2)
library(plotly)
library(pheatmap)

server <- function(input, output, session) {
  
  classic_purple_gray <- c(
    "#7C75C4", "#A48FBD", "#C9A4BA", "#E3B4C2",
    "#EBD4D0", "#7E7C48", "#949773", "#A39F53",
    "#BFA7D5", "#D0B7E1", "#CDBDBC", "#DDD2C7"
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
    withProgress(message = "Loading Heatmap...", value = 0.5, {
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
        color = colorRampPalette(c("blue", "white", "red"))(50),
        main = paste0("Top ", n_genes, " Variable Genes (50 Random Tumor Samples)"),
        show_colnames = FALSE,
        show_rownames = TRUE,
        fontsize_row = 10
      )
    })
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
            type = 'pie',
            textposition = 'inside',
            textinfo = 'label+percent',
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
  
  # Download Data
  merged_data_reactive <- reactive({
    req(input$cancer_type) 
    
    obj <- get_cancer_object(input$cancer_type)
    
    return(obj$merged_data)
  })
  
  output$dl_data_preview <- renderTable({
    df <- merged_data_reactive()
    head(df, 20)
  })
  
  output$download_data_btn <- downloadHandler(
    filename = function() {
      paste0(input$cancer_type, "_Merged_Data.csv")
    },
    content = function(file) {
      write.csv(merged_data_reactive(), file, row.names = FALSE)
    }
  )
}

shinyServer(server)