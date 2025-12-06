library(shiny)
library(shinycssloaders)
library(plotly)

# DE Gene Sidebar
deGeneSidebar <- function() {
  sidebarPanel(
    width = 4,
    div(
      class = "panel-with-header",
      div(
        class = "panel-header",
        h4("DE Gene Inputs")
      ),
      div(
        class = "panel-body",
        # Subtitle section
        tagList(
          tags$small("Configure parameters for differential expression benchmarking."),
          tags$hr()
        ),
        
        # Inputs
        tags$div(
          class = "sidebar-section",
          tags$h4("Tumor site"),
          selectInput(
            inputId = "cancer_type_de",
            label = "Cancer type",
            choices = c("Thyroid", "Liver", "Kidney", "Breast", "Colon", "Bronchus and Lung"),
            selected = "Thyroid"
          )
        ),
        
        tags$hr(),
        
        tags$div(
          class = "sidebar-section",
          tags$h4("Sampling"),
          sliderInput(
            inputId = "de_sample_size",
            label = "Sample size",
            min = 5,
            max = 25,
            value = 10,
            step = 1
          )
        ),
        
        tags$hr(),
        
        tags$div(
          class = "sidebar-section",
          tags$h4("Signal strength"),
          sliderInput(
            inputId = "de_pct",
            label = "Percentage of DE genes (%)",
            min = 5,
            max = 15,
            value = 10,
            step = 1
          )
        ),
        
        tags$hr(),
        
        tags$div(
          class = "sidebar-section",
          tags$h4("Dispersion Spread ϕ"),
          sliderInput(
            inputId = "dispersion",
            label = "Dispersion (phi ϕ)",
            min = 0.1,
            max = 0.5,
            value = 0.3,
            step = 0.1
          )
        ),
        
        div(
          class = "sidebar-actions",
          actionButton(
            inputId = "run_de",
            label = "Run DE Genes Detection",
            width = "100%"
          )
        )
      )
    )
  )
}

selected_genes <- c(
  "COL11A1", "MMP1", "MMP7", "LCN2", "LGR5",
  "PROM1", "IGF2", "EREG", "CXCL13", "VTCN1",
  "LCN2", "ADH1B", "PCK1", "AKR1B10", "CA4",
  "ROS1", "SFRP1", "TFF1","XIST")

# Survival Sidebar
survivalSidebar <- function() {
  sidebarPanel(
    width = 4,
    div(
      class = "panel-with-header",
      div(
        class = "panel-header",
        h4("Survival Inputs")
      ),
      div(
        class = "panel-body",
        # Subtitle section
        tagList(
          tags$small("Configure demographic and genomic parameters for the survival analysis."),
          tags$hr()
        ),
        
        # Inputs
        tags$div(
          class = "sidebar-section",
          tags$h4("Tumor site"),
          selectInput(
            inputId = "cancer_type_surv",
            label   = "Cancer type",
            choices = c("Thyroid", "Liver", "Kidney", "Breast", "Colon", "Bronchus and Lung"),
            selected = "Bronchus and Lung"
          )
        ),
        
        tags$hr(),
        
        tags$div(
          class = "sidebar-section",
          tags$h4("Demographics"),
          selectInput(
            inputId = "gender_input",
            label = "Sex",
            choices = c("Female", "Male"),
            selected = "Female"
          ),
          sliderInput(
            inputId = "age_input",
            label = "Age",
            min = 18,
            max = 99,
            value = 55
          )
        ),
        
        tags$hr(),
        
        tags$div(
          class = "sidebar-section",
          tags$h4("Genomics"),
          selectizeInput(
            inputId  = "selected_gene",
            label = "Gene",
            choices = selected_genes,
            selected = "COL11A1",
            multiple = FALSE
          ),
          sliderInput(
            inputId = "expr_treated",
            label = "Gene Expression Level (Log 2)",
            min = 0,
            max = 12,
            value = 6,
            step = 0.5
          )
        ),
        
        div(
          class = "sidebar-actions",
          actionButton(
            inputId = "run_sim",
            label = "Run Survival Analysis",
            width = "100%"
          )
        )
      )
    )
  )
}

# Download Sidebar 
downloadSidebar <- function() {
  sidebarPanel(
    width = 4,
    div(
      class = "panel-with-header",
      div(
        class = "panel-header",
        h4("Data Set Filters")
      ),
      div(
        class = "panel-body",
        # Subtitle section
        tagList(
          tags$small("Select a cancer type to preview and download the merged and cleaned dataset."),
          tags$hr()
        ),
        
        # Inputs
        tags$div(
          class = "sidebar-section",
          tags$h4("Tumor site"),
          selectInput(
            inputId = "cancer_type_dl",
            label = "Cancer type",
            choices = c("Thyroid", "Liver", "Kidney", "Breast", "Colon", "Bronchus and Lung"),
            selected = "Thyroid"
          )
        )
      )
    )
  )
}

ui <- tagList(
  # Global CSS
  tags$head(
    tags$style(HTML("
      /* Spacing */
      .container-fluid {
        max-width: 1500px;
      }
      
      /* NavBar */
      .navbar.navbar-default {
        background-color: #7A658A;
        border-color: #89709E;
      }
      
      .navbar-default .navbar-brand,
      .navbar-default .navbar-nav > li > a {
        color: #ffffff;
      }
      
      .navbar-default .navbar-brand:hover,
      .navbar-default .navbar-nav > li > a:hover,
      .navbar-default .navbar-nav > li > a:focus {
        color: #ffcc66;
      }
      
      .nav-tabs > li > a {
      color: #5B469E;
      font-weight: 500;
    }
      /* body */
      body {
        font-family: 'Verdana', Arial, sans-serif;
        font-size: 14px;
        padding-bottom: 70px;
      }

      /* human body container */
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

      /* Panel */
      .panel-with-header {
        border: 1px solid #e3e3e3;
        border-radius: 4px;
        margin-bottom: 15px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
      }

      .panel-header {
        background-color: #7A658A;
        color: white;
        padding: 12px 15px;
        margin: 0;
        border-radius: 4px 4px 0 0;
      }

      .panel-header h4 {
        margin: 0;
        font-size: 16px;
        font-weight: 600;
        color: white;
      }

      .panel-body {
        background-color: #f9f9f9;
        padding: 15px;
      }

      /* Sidebar styling */
      .sidebar-panel {
        padding: 0;
      }
      
      .sidebar-section h4 {
        font-size: 14px;
        font-weight: 600;
        margin-top: 6px;
        margin-bottom: 6px;
        color: #333;
      }
      
      .sidebar-section label,
      .panel-body label {
        font-size: 12px;
        font-weight: 550;
      }
      
      .panel-body hr {
        margin: 0px 0;
        border-top: 1px solid #7A658A;
      }
      
      .panel-body small {
        color: #666;
        display: block;
        margin-bottom: 8px;
      }
      
      /* Sidebar sections */
      .sidebar-section {
        padding-top: 4px;
        margin-bottom: 6px;
      }
      
      /* Action buttons */
      .sidebar-actions {
        margin-top: 12px;
        padding-top: 12px;
        border-top: 1px solid #7A658A;
      }
      
      .sidebar-actions .btn,
      .sidebar-section .btn {
        width: 100%;
        font-weight: 400;
        margin-top: 6px;
        border-radius: 3px;
      }
      
      /* Slider Fill */
      .irs--shiny .irs-bar,
      .irs--shiny .irs-handle,
      .irs--shiny .irs-from,
      .irs--shiny .irs-to,
      .irs--shiny .irs-single {
        background: #7A658A;
      }
      
      .irs--shiny .irs-bar {
        background: #7A658A;
        border-top: 1px solid #7A658A;
        border-bottom: 1px solid #7A658A;
      }
      
      .irs--shiny .irs-from:before,
      .irs--shiny .irs-to:before,
      .irs--shiny .irs-single:before {
        border-top-color: #7A658A;
      }
      
      /* Button styling */
      .sidebar-actions .btn-default,
      .sidebar-section .btn-default {
        background-color: #7A658A;
        border-color: #7A658A;
        color: #ffffff;
      }
      
      .sidebar-actions .btn-default:hover,
      .sidebar-actions .btn-default:focus,
      .sidebar-section .btn-default:hover,
      .sidebar-section .btn-default:focus {
        background-color: #614d6e;
        border-color: #614d6e;
        color: #ffffff;
      }
      
      .well {
        background-color: transparent;
        border: none;
        box-shadow: none;
        padding: 0px; 
      }
      
      .app-footer {
        background-color: #7A658A;
        color: white;
        padding: 3px;
        text-align: center;
        margin-top: 0;
        width: 100%;
        position: fixed;
        bottom: 0;
        left: 0;
      }
    "))
  ),
  
  navbarPage(
    title = "The Pan-Cancer Demographic-Genomic Simulator",
    
    # Body Map Page
    tabPanel(
      title = "🩻 Body Map",
      fluidPage(
        h3("Body Map"),
        p("Hover over and click different regions to explore different cancers in their demographics and TOP20 gene heatmap."),
        
        fluidRow(
          # Body Map
          column(
            width = 4,
            div(
              id = "body-wrapper",
              tags$img(
                id = "body-image",
                src = "human_body.svg",
                alt = "Human body diagram",
                class = "img-responsive"
              ),
              
              # SVG overlay
              tags$div(
                id = "body-svg-overlay",
                HTML('
              <svg viewBox="0 0 226 917" xmlns="http://www.w3.org/2000/svg">
                <!-- Thyroid -->
                <ellipse class="organ-path" data-organ="Thyroid" id="organ-thyroid"
                  cx="113" cy="150" rx="15" ry="10" 
                  fill="#9ECAE1" stroke="#7AA3C1" stroke-width="2"
                  onclick="Shiny.setInputValue(\'selected_organ\', \'Thyroid\', {priority: \'event\'})" />
                  
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
                  fill="#BCBD22" stroke="#9A9B1C" stroke-width="2"
                  onclick="Shiny.setInputValue(\'selected_organ\', \'Bronchus and Lung\', {priority: \'event\'})" />

                <!-- Breast -->
                <ellipse class="organ-path" data-organ="Breast" id="organ-breast"
                  cx="150" cy="240" rx="30" ry="22" 
                  fill="#DD3497" stroke="#B32878" stroke-width="2"
                  onclick="Shiny.setInputValue(\'selected_organ\', \'Breast\', {priority: \'event\'})" />
                
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
                  fill="#E7969C" stroke="#C7767C" stroke-width="2"
                  onclick="Shiny.setInputValue(\'selected_organ\', \'Liver\', {priority: \'event\'})" />
                
               <!-- Kidney -->
                <ellipse class="organ-path" data-organ="Kidney" id="organ-kidney-left"
                  cx="80" cy="340" rx="12" ry="18" 
                  fill="#CE6DBD" stroke="#AE4D9D" stroke-width="2"
                  onclick="Shiny.setInputValue(\'selected_organ\', \'Kidney\', {priority: \'event\'})" />
                <ellipse class="organ-path" data-organ="Kidney" id="organ-kidney-right"
                  cx="146" cy="340" rx="12" ry="18" 
                  fill="#CE6DBD" stroke="#AE4D9D" stroke-width="2"
                  onclick="Shiny.setInputValue(\'selected_organ\', \'Kidney\', {priority: \'event\'})" />
                
                <!-- Colon -->
                <path class="organ-path" data-organ="Colon" id="organ-colon"
                  d="M 85,370 
                  L 85,430 
                  Q 85,450 113,450 
                  Q 141,450 141,430 
                  L 141,370 
                  Q 141,360 113,360 
                  Q 85,360 85,370 
                  Z"
                  fill="#17BECF" stroke="#129EAF" stroke-width="2"
                  onclick="Shiny.setInputValue(\'selected_organ\', \'Colon\', {priority: \'event\'})" />
              </svg>
            ')
              )
            )
          ),
          
          # body map plots and stats
          column(
            width = 8,
            
            # when nothing is selected 
            conditionalPanel(
              condition = "input.selected_organ == null || input.selected_organ == ''",
              div(
                class = "panel-with-header",
                div(class = "panel-header", h4("Body Map Explorer")),
                div(
                  class = "panel-body",
                  p(em("Click on a highlighted region in the body to explore more!"))
                )
              )
            ),
            
            # when an organ is selected
            conditionalPanel(
              condition = "input.selected_organ != null && input.selected_organ != ''",
              
              fluidRow(
                column(
                  width = 12,
                  div(
                    class = "panel-with-header",
                    div(class = "panel-header", 
                        h4(textOutput("organ_title", inline = TRUE))),
                    div(class = "panel-body", 
                        withSpinner(
                          uiOutput("organ_intro"),
                          type = 4,
                          color = "#7A658A",
                          size = 0.7,
                          proxy.height = "100px"
                          ))
                  )
                )
              ),
              
              fluidRow(
                column(
                  width = 12,
                  div(
                    class = "panel-with-header",
                    div(class = "panel-header", h4("Analysis Plots")),
                    div(
                      class = "panel-body",
                      tabsetPanel(
                        id = "analysis_tabs",
                        
                        # Tab 1: Heatmap
                        tabPanel(
                          title = "Top 20 Genes Expression",
                          br(),
                          withSpinner(plotOutput("plot_heatmap", height = "500px"),
                                      type = 4,
                                      color = "#7A658A",
                                      size = 1)
                        ),
                        
                        # Tab 2: Race
                        tabPanel(
                          title = "Race Breakdown",
                          br(),
                          plotlyOutput("plot_race", height = "500px")
                        ),
                        
                        # Tab 3: Sex
                        tabPanel(
                          title = "Sex Breakdown",
                          br(),
                          plotlyOutput("plot_gender", height = "500px")
                        ),
                        
                        # Tab 4: Age
                        tabPanel(
                          title = "Age Breakdown",
                          br(),
                          plotlyOutput("plot_age", width = "800px", height = "500px")
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    ),
    
    # DE Gene Analysis
    tabPanel(
      title = "🧬 DE Genes Detection",
      sidebarLayout(
        deGeneSidebar(),
        mainPanel = mainPanel(
          width = 8,
          h3("DE Genes Detection Benchmark"),
          plotOutput("de_genes_detection")
        )
      )
    ),
    
    # Survival Analysis
    tabPanel(
      title = "📉 Survival Analysis",
      sidebarLayout(
        survivalSidebar(),
        mainPanel = mainPanel(
          width = 8,
          tabsetPanel(
            id = "surv_tabs",
            tabPanel(
              title = "Survival Curves",
              uiOutput("surv_curv_ui")
            ),
            tabPanel(
              title = "Survival Gains",
              uiOutput("surv_gain_ui")
            )
          ),
          plotlyOutput("surv_plot", height = "450px")
        )
      )
    )
    ,
    
    # Data Set
    tabPanel(
      title = "🗃️ Data Set",
      sidebarLayout(
        downloadSidebar(),
        
        mainPanel = mainPanel(
          width = 8,
          h3("Download Merged Data"),
          p("Download the fully merged clinical and genomic dataset for the cancer type selected in the sidebar."),
          br(),
          
          downloadButton("download_data_btn", "Merged Data (.csv)", class = "btn-default"),
          
          tags$hr(),
          
          h4("Preview (Top 20 Rows)"),
          div(
            style = "overflow-x: scroll; border: 1px solid #eee; padding: 10px;", 
            withSpinner(
              tableOutput("dl_data_preview"), 
              type = 4,
              color = "#7A658A",
              size = 1
            )
          )
        )
      )
    ),
    
    footer = tags$div(
      class = "app-footer",
      p("© 2025 Pan-Cancer Simulator | Developed with R Shiny & Plotly | Liuhan Ke, Ziyi Ou, Hailin Zhang")
    )
  )
)

shinyUI(ui)