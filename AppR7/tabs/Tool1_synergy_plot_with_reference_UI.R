####UI-2.1  Synergy Plot (Reference)  ####
tabPanel("Synergy Plot",
         sidebarLayout(fluid = FALSE, position="left",sidebarPanel(
           br(),
           tags$div(
             p("1. Disease Signature", style = "text-align: center;margin:0px;color: #34495e;font-weight: bold;font-size: 20px;"),
             hr(style="border-top-style: solid;border-top-width:2px;margin:12px;"),
             selectInput(inputId = "disease",
                          label = "Select a Disease Signature:",
                          choices = c("TCGA Bladder Urothelial Carcinoma (BLCA)", 
                                      "TCGA Breast Invasive Carcinoma (BRCA)", 
                                      "TCGA Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESC)", 
                                      "TCGA Cholangiocarcinoma (CHOL)", 
                                      "TCGA Colon Adenocarcinoma (COAD)", 
                                      "TCGA Esophageal Carcinoma (ESCA)", 
                                      "TCGA Glioblastoma Multiforme (GBM)", 
                                      "TCGA Head and Neck Squamous Cell Carcinoma (HNSC)", 
                                      "TCGA Kidney Chromophobe (KICH)", 
                                      "TCGA Kidney Renal Clear Cell Carcinoma (KIRC)", 
                                      "TCGA Kidney Renal Papillary Cell Carcinoma (KIRP)", 
                                      "TCGA Liver Hepatocellular Carcinoma (LIHC)", 
                                      "TCGA Lung Adenocarcinoma (LUAD)", 
                                      "TCGA Lung Squamous Cell Carcinoma (LUSC)", 
                                      "TCGA Pheochromocytoma and Paraganglioma (PCPG)", 
                                      "TCGA Prostate Adenocarcinoma (PRAD)", 
                                      "TCGA Rectum Adenocarcinoma (READ)", 
                                      "TCGA Stomach Adenocarcinoma (STAD)", 
                                      "TCGA Thyroid Carcinoma (THCA)", 
                                      "TCGA Uterine Corpus Endometrial Carcinoma (UCEC)", 
                                      "TCGA Glioblastoma signature (GBM.sig)", 
                                      "TCGA Colon signature (CRC.sig)", 
                                      "TCGA Breast signature (BRCA.sig)", 
                                      "PDX GBM Group 1 (PDX.GBM1)", 
                                      "PDX GBM Group 2 (PDX.GBM2)", 
                                      "PDX GBM Group 3 (PDX.GBM3)", 
                                      "PDX GBM Group 4 (PDX.GBM4)")
                         ),
             helpText("OR",style="text-align: center;margin-top:0px;margin-bottom: 6px;"),
             textAreaInput("Ref_Sig_Text_Box", "Paste a Disease Signature", placeholder="Gene, Expression Value"),
             style="border-width: 2px;padding:15px;border-radius: 5px; border-color: #34495e;background-color: white;border-style: solid;"),
           
           br(),
           br(),
           
           tags$div(
             p("2. Reference Drug", style = "text-align: center;margin:0px;color: #34495e;font-weight: bold;font-size: 20px;"),
             hr(style="border-top-style: solid;border-top-width:2px;margin:12px;"),
             uiOutput("ui",style="margin-bottom: 0px;"),
             helpText("OR",style="text-align: center;margin-top:0px;margin-bottom: 6px;"),
             textAreaInput("Ref_Sig_Text_Box", "Paste a Reference Signature", placeholder="Gene, Expression Value"),
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
           helpText("Filter out the lowest 'n' percentile of the gene consensus scores from the Reference Drug Signature"),
           
           br(),
           selectInput(inputId = "bioactivities",
                       label = "Choose a measure of biological activity:",
                       choices = c("IC50", "Kd", "Ki", "Potency"), selected = "IC50"),
           uiOutput("ui_bioact"),
           helpText("Choose a molecular drug target of interest to compare drug activity."),
           br() 
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
