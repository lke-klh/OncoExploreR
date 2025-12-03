library(shiny)

simulationSidebar <- function() {
  sidebarPanel(
    h4("Simulation Inputs"),
    
    selectInput(
      inputId = "gender",
      label   = "Gender",
      choices = c("Female", "Male", "Non-binary", "Other"),
      selected = "Female"
    ),
    
    selectInput(
      inputId = "race",
      label   = "Race",
      choices = c("White", "Black", "Asian", "Hispanic", "Other"),
      selected = "White"
    ),
    
    sliderInput(
      inputId = "age",
      label   = "Age",
      min     = 0,
      max     = 100,
      value   = 50
    ),
    
    selectizeInput(
      inputId  = "genes",
      label    = "Select gene(s)",
      choices  = c("TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "PIK3CA"),
      selected = "TP53",
      multiple = TRUE
    ),
    
    sliderInput(
      inputId = "gene_expr",
      label   = "Gene expression level",
      min     = 0,
      max     = 10,
      value   = 5,
      step    = 0.1
    ),
    
    selectInput(
      inputId = "cancer_type",
      label   = "Cancer type",
      choices = c("Thyroid", "Liver", "Kidney", "Breast", "Colon", "Bronchus and Lung"),
      selected = "Thyroid"
    ),
    
    br(),
    actionButton(
      inputId = "run_sim",
      label   = "Run Simulation"
    )
  )
}

ui <- tagList(
  # Global CSS (applies to all pages)
  tags$head(
    tags$style(HTML("
      /* Spacing */
      .container-fluid {
        max-width: 1200px;
      }

      /* Body container with image and hover overlays */
      #body-wrapper {
        position: relative;
        max-width: 110px;
        margin: 0 auto 10px auto;
      }

      #body-wrapper img {
        width: 100%;
        height: auto;
        display: block;
      }
      
      #body-svg-overlay {
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        pointer-events: none;
      }
      
      #body-svg-overlay svg {
        width: 100%;
        height: 100%;
      }

      .organ-path {
        opacity: 0.4;
        cursor: pointer;
        transition: opacity 0.3s ease, filter 0.3s ease;
        pointer-events: auto;
      }

      .organ-path:hover {
        opacity: 0.85;
        filter: brightness(1.3) saturate(1.5);
      }
      
      .organ-path.active {
        opacity: 0.85;
        filter: brightness(1.3) saturate(1.5);
      }
      
      /* Tooltip styling */
      .body-tooltip {
        position: absolute;
        opacity: 0;
        background-color: white;
        padding: 10px;
        box-shadow: 0 2px 5px 0 rgba(0,0,0,0.16), 0 2px 10px 0 rgba(0,0,0,0.12);
        border-radius: 5px;
        border: 1px solid rgba(40, 40, 40);
        pointer-events: none;
        z-index: 1000;
        transition: opacity 0.2s;
      }
      
      .body-tooltip .organ-name {
        color: #bb0e3d;
        font-weight: bold;
        margin-bottom: 5px;
      }
      
      .body-tooltip .organ-info {
        font-size: 12px;
        color: rgb(20, 20, 20);
      }
    "))
  ),
  
  navbarPage(
    title = "The Pan-Cancer Demographic-Genomic Simulator",
    
    # ---- Body Map page (no sidebar) ----
    tabPanel(
      title = "Body Map",
      fluidPage(
        h3("Body Map"),
        p("Hover over different regions of the body to highlight them."),
        
        div(
          id = "body-wrapper",
          tags$img(
            id   = "body-image",
            src  = "human_body.svg",
            alt  = "Human body diagram",
            class = "img-responsive"
          ),
          
          # SVG overlay
          tags$div(
            id = "body-svg-overlay",
            HTML('
              <svg viewBox="0 0 226 917" xmlns="http://www.w3.org/2000/svg">
                <!-- Thyroid (neck area - small oval) -->
                <ellipse class="organ-path" data-organ="Thyroid" id="organ-thyroid"
                  cx="113" cy="150" rx="15" ry="10" 
                  fill="#9ECAE1" stroke="#7AA3C1" stroke-width="2"/>
                  
                <!-- Bronchus and Lung -->
                <path class="organ-path" data-organ="Bronchus and Lung" id="organ-lung"
                  d="M 105,245
                     C 103,236 102,228 102,220
                     C 102,207 100,190 88,190
                     C 72,190 63,220 63,250
                     C 63,252 63,270 66,274
                     C 68,277 71,278 73,278
                     C 79,278 82,276 86,272
                     C 92,267 99,261 105,257
                     C 110,254 108,251 105,245 Z"
                  fill="#BCBD22" stroke="#9A9B1C" stroke-width="2"/>

                <!-- Breast -->
                <ellipse class="organ-path" data-organ="Breast" id="organ-breast"
                  cx="150" cy="240" rx="30" ry="22" 
                  fill="#DD3497" stroke="#B32878" stroke-width="2"/>
                
                <!-- Liver -->
                <path class="organ-path" data-organ="Liver" id="organ-liver"
                  d="M 83,305 
                     C 80,302 79,299 80,296
                     C 81,293 85,291 90,290
                     C 97,289 108,288 119,289
                     C 128,290 135,291 140,293
                     C 145,295 148,297 149,300
                     C 150,303 148,306 145,309
                     C 141,313 134,316 125,319
                     C 114,322 101,323 91,323
                     C 84,323 79,322 76,319
                     C 74,317 73,314 74,311
                     L 83,305 Z"
                  fill="#E7969C" stroke="#C7767C" stroke-width="2"/>
                
               <!-- Kidney -->
                <ellipse class="organ-path" data-organ="Kidney" id="organ-kidney-left"
                  cx="80" cy="340" rx="12" ry="18" 
                  fill="#CE6DBD" stroke="#AE4D9D" stroke-width="2"/>
                <ellipse class="organ-path" data-organ="Kidney" id="organ-kidney-right"
                  cx="146" cy="340" rx="12" ry="18" 
                  fill="#CE6DBD" stroke="#AE4D9D" stroke-width="2"/>
                
                <!-- Colon -->
                <path class="organ-path" data-organ="Colon" id="organ-colon"
                  d="M 85,370 L 85,430 Q 85,450 113,450 Q 141,450 141,430 L 141,370 Q 141,360 113,360 Q 85,360 85,370 Z"
                  fill="#17BECF" stroke="#129EAF" stroke-width="2"/>
              </svg>
            ')
          )
        ),
        br(),
        conditionalPanel(
          condition = "input.selected_organ != null && input.selected_organ != ''",
          fluidRow(
            column(
              width = 4,
              wellPanel(
                h4(textOutput("organ_title")),
                uiOutput("organ_intro")
              )
            ),
            
            column(
              width = 8,
              wellPanel(
                h4("Selected Cancer: Plots & Summary"),
                plotOutput("organ_plot"),
                br(),
                h5("Summary statistics"),
                verbatimTextOutput("organ_stats")
              )
            )
          )
        ),
      )
    ),
    
    # Survival Analysis
    tabPanel(
      title = "Survival Analysis",
      sidebarLayout(
        sidebarPanel = simulationSidebar(),
        mainPanel(
          tabsetPanel(
            id = "surv_tabs",
            
            tabPanel(
              title = "Survival Analysis",
              h3("Survival Analysis Plot"),
              plotOutput("surv_plot")
            ),
            
            tabPanel(
              title = "Summary Stats",
              h3("Summary Statistics"),
              verbatimTextOutput("summary_stats")
            )
          )
        )
      )
    ),
    
    # Distribution Page
    tabPanel(
      title = "Distribution",
      sidebarLayout(
        sidebarPanel = simulationSidebar(),
        mainPanel(
          h3("Distribution of Simulated Data"),
          plotOutput("dist_plot")
        )
      )
    )
  )
)

shinyUI(ui)
