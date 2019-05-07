####UI-2.1  Synergy Plot (Reference)  ####
tabPanel("Synergy Plot",
         sidebarLayout(fluid = FALSE, position="left",sidebarPanel(
           br(),
           tags$div(
             p("1. Disease Signature", style = "text-align: center;margin:0px;color: #34495e;font-weight: bold;font-size: 20px;"),
             hr(style="border-top-style: solid;border-top-width:2px;margin:12px;"),
             selectInput(inputId = "T1_Disease",
                         label = "Select a Disease Signature:",
                         choices = c("Glioblastoma TCGA (GBM)","Colon TCGA (CRC)", "Breast TCGA (BRCA)",
                                     "PDX GBM Group 1", "PDX GBM Group 2", "PDX GBM Group 3", "PDX GBM Group 4")),
             helpText("OR",style="text-align: center;margin-top:0px;margin-bottom: 6px;"),
             fluidRow(
               column(8,  h5("Paste a Disease Signature",style="font-weight: bold") )), 
             fluidRow(
               
               column(8, textAreaInput("T1_Disease_Text_Box",label=NULL, placeholder="Gene, Expression Value")),
               column(2, actionButton("T1_go_disease_signature", "Go"))),
             column(1),
             style="border-width: 2px;padding:15px;border-radius: 5px; border-color: #34495e;background-color: white;border-style: solid;"),
           
           br(),
           br(),
           
           tags$div(
             p("2. Reference Drug", style = "text-align: center;margin:0px;color: #34495e;font-weight: bold;font-size: 20px;"),
             hr(style="border-top-style: solid;border-top-width:2px;margin:12px;"),
             uiOutput("ui",style="margin-bottom: 0px;"),
             helpText("OR",style="text-align: center;margin-top:0px;margin-bottom: 6px;"),
             fluidRow(
               column(8,  h5("Paste a Reference Signature",style="font-weight: bold") )), 
             fluidRow(
               
               column(8, textAreaInput("T1_Reference_Text_Box",label=NULL, placeholder="Gene, Expression Value")),
               column(2, actionButton("T1_go_reference_signature", "Go"))),
             column(1),
             style="border-width: 2px;padding:15px;border-radius: 5px; border-color: #34495e;background-color: white;border-style: solid;"),
           
           
           br(),
           br(),
           hr(style="margin:5px;border-top-style: dotted;border-top-width:2px;border-top-color:#34495e;"),
           p("Options:", style = "text-align: center;margin:0px;font-weight: bold;font-size: 20px;"),
           hr(style="margin:5px;border-top-style: dotted;border-top-width:2px;border-top-color:#34495e;"),
           br(),
           selectInput(inputId = "L1000_Dataset",
                       label = "Choose an L1000 Dataset:",
                       choices = c("LINCS L1000 Dec 2015", "LINCS L1000 March 2017"),selected = "LINCS L1000 Dec 2015"),
           
           sliderInput(inputId = "bins",
                       label = "",
                       min = 1,
                       max = 100,
                       value = 33),
           helpText("Filter out the lowest 'n' percentile of the gene consensus scores from the Reference Drug Signature")
         ),
         
         
         mainPanel(fluid = FALSE,
                   
                   tags$ul(
                     #tags$li(textOutput("selected_var")),
                     tags$li(textOutput("selected_var2")),
                     tags$li(textOutput("selected_var3"))
                   ),
                   plotlyOutput("plot1"),
                   br(),
                   br(),
                   br(),
                   br(),
                   downloadButton("downloadData", "Download Synergy Plot Data"),
                   downloadButton("downloadData_Drug_Sig", "Download Drug Signatures"),
                   DT::dataTableOutput("view1")
                   #tableOutput("table1")
         ))
         
)
