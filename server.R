library(shiny)
library(tidyverse)
library(ggplot2)
library(plotly)

server <- function(input, output, session) {
  
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
  
  output$organ_plot <- renderPlot({
    req(input$selected_organ)
    plot(1:10, rnorm(10), main = paste("Placeholder plot for", input$selected_organ))
  })
}



shinyServer(server)
