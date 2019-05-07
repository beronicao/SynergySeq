####UI-3.1  Single Cell Synergy ####
tabPanel("Single Cell RNAseq",
         sidebarLayout(fluid = FALSE, position="left",sidebarPanel(
           br(),
           tags$div(
             p("1. Load Single Cell Dataset", style = "text-align: center;margin:0px;color: #34495e;font-weight: bold;font-size: 20px;"),
             hr(style="border-top-style: solid;border-top-width:2px;margin:12px;"),
             selectInput(inputId = "T3_Single_Cell_File",
                         label = "Choose a Single Cell RNAseq Dataset",
                         choices = c("Glioblastoma")),
             helpText("Or"),
             fileInput("T3_Single_Cell_File2", "Upload a Single Cell RNAseq Dataset",
                       multiple = FALSE,
                       accept = c(".txt")),
             downloadButton("T3_Single_Cell_File_Example", "Download Template"),

             style="border-width: 2px;padding:15px;border-radius: 5px; border-color: #34495e;background-color: white;border-style: solid;")
           
           # br(),
           # 
           # tags$div(
           #   p("2. Choose Compound", style = "text-align: center;margin:0px;color: #34495e;font-weight: bold;font-size: 20px;"),
           #   hr(style="border-top-style: solid;border-top-width:2px;margin:12px;"),
           #   #uiOutput("T3_ui"),
           #   style="border-width: 2px;padding:15px;border-radius: 5px; border-color: #34495e;background-color: white;border-style: solid;")
           # 
           
           
         ),
         mainPanel(fluid = FALSE,
                   br(),
                   br(),
                   uiOutput("T3_ui"),
                   plotlyOutput("T3_plot1"))
         
         
         )
         
         
         
)
