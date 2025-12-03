library(shiny)
library(tidyverse)
library(ggplot2)
library(plotly)

server <- function(input, output, session) {
  
  output$organ_title <- renderText({
    req(input$selected_organ)
    paste("Selected cancer site:", input$selected_organ)
  })
  
  # Intro text on the left
  output$organ_intro <- renderUI({
    req(input$selected_organ)
    
    text <- switch(
      input$selected_organ,
      "Thyroid"            = "Thyroid cancer arises from the thyroid gland...",
      "Bronchus and Lung"  = "Lung cancer is commonly divided into small cell and non-small cell types...",
      "Breast"             = "Breast cancer is one of the most common malignancies among women...",
      "Liver"              = "Primary liver cancer often refers to hepatocellular carcinoma...",
      "Kidney"             = "Kidney cancer includes clear cell, papillary, and other histologic subtypes...",
      "Colon"              = "Colorectal cancer typically develops from adenomatous polyps...",
      "Selected site."
    )
    
    tagList(
      p(text),
      p(em("You can replace this text with any detailed description later."))
    )
  })
  
  output$organ_plot <- renderPlot({
    req(input$selected_organ)
    plot(1:10, rnorm(10), main = paste("Placeholder plot for", input$selected_organ))
  })
  
  # Placeholder summary stats
  output$organ_stats <- renderPrint({
    req(input$selected_organ)
    list(
      cancer_site = input$selected_organ,
      n_patients  = 1234,
      median_age  = 65,
      notes       = "Replace this list with real summary stats from your backend."
    )
  })
}



shinyServer(server)
