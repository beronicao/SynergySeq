tabPanel("Synergy Plot (No Reference Drug)",
         sidebarLayout(fluid = FALSE, position="left",sidebarPanel(
           br(),
           tags$div(
             p("1. Disease Signature", style = "text-align: center;margin:0px;color: #34495e;font-weight: bold;font-size: 20px;"),
             hr(style="border-top-style: solid;border-top-width:2px;margin:12px;"),
             selectInput(inputId = "T2_Disease",
                         label = "Select a Disease Signature:",
                         choices = c("","Glioblastoma TCGA (GBM)","Colon TCGA (CRC)", "Breast TCGA (BRCA)",
                                     "PDX GBM Group 1", "PDX GBM Group 2", "PDX GBM Group 3", "PDX GBM Group 4")),
             style="border-width: 2px;padding:15px;border-radius: 5px; border-color: #34495e;background-color: white;border-style: solid;"),
           br(),
           br(),
           hr(style="margin:5px;border-top-style: dotted;border-top-width:2px;border-top-color:#34495e;"),
           p("Options:", style = "text-align: center;margin:0px;font-weight: bold;font-size: 20px;"),
           hr(style="margin:5px;border-top-style: dotted;border-top-width:2px;border-top-color:#34495e;"),
           br(),
           selectInput(inputId = "T2_Simil_metric",
                       label = "Select a Drug Similarity Metric:",
                       choices = c("Pearson Correlation"),selected = "Pearson Correlation"),
                       #choices = c("Similarity Ratio","Weighted Similarity Ratio","Pearson Correlation"),selected = "Pearson Correlation"),
           br()
           
         ),
         
         mainPanel(fluid = FALSE,
                   br(),
                   
                   
                   textOutput("T2_selected_var2"),
                   plotlyOutput("T2_plot1"),
                   br(),
                   tableOutput("T2_value"),
                   br(),
                   br(),
                   br(),
                   DT::dataTableOutput("T2_view1")
         )
         
         ))