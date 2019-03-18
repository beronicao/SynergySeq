#App7



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### 1. LOADING ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
library(ggplot2)
library(shiny)
library(plotly)
library(DT)
library(visNetwork)
library(reshape2)
library(shinythemes)
library(RColorBrewer)

SM_MOA <- read.csv(file="data/L1000_SM_MOA.csv", header=TRUE)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### 2. UI ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
ui <- navbarPage(  theme = shinytheme("flatly"),responsive = TRUE,
                   title = div(
                              img(src = "Logos_3.png",height="30",align="left",position="fixed",hspace="10"),
                              "SynergySeq"
                              ),
                  
                   ### UI Introduction ####
                   #source("tabs/introduction_UI.R",  local = TRUE)$value,
                  
                   
                   tabPanel("Tools", tabsetPanel(
                     
                            ### UI Tool1 ####
                            source("tabs/Tool1_synergy_plot_with_reference_UI.R",  local = TRUE)$value,
                            
                            
                            ### UI Tool2 ####
                            tabPanel("Synergy Plot (No Reference Drug)",
                                     sidebarLayout(fluid = FALSE, position="left",sidebarPanel(
                                                  br(),
                                                  tags$div(
                                                          p("1. Disease Signature", style = "text-align: center;margin:0px;color: #34495e;font-weight: bold;font-size: 20px;"),
                                                          hr(style="border-top-style: solid;border-top-width:2px;margin:12px;"),
                                                          selectInput(inputId = "T2_Disease",
                                                            label = "Select a Disease Signature:",
                                                            choices = c("",
                                                                      "TCGA Bladder Urothelial Carcinoma (BLCA)", 
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
                                                          fluidRow(
                                                            column(8,  h5("Paste a Disease Signature",style="font-weight: bold") )), 
                                                          fluidRow(
                                                            
                                                            column(8, textAreaInput("T2_Disease_Text_Box",label=NULL, placeholder="Gene, Expression Value")),
                                                            column(2, actionButton("T2_go_disease_signature", "Go"))),
                                                          column(1),
                                                          style="border-width: 2px;padding:15px;border-radius: 5px; border-color: #34495e;background-color: white;border-style: solid;"),
                                                    br(),
                                                    br(),
                                                    hr(style="margin:5px;border-top-style: dotted;border-top-width:2px;border-top-color:#34495e;"),
                                                    p("Options:", style = "text-align: center;margin:0px;font-weight: bold;font-size: 20px;"),
                                                    hr(style="margin:5px;border-top-style: dotted;border-top-width:2px;border-top-color:#34495e;"),
                                                    br(),
                                                    sliderInput(inputId = "T2_bins",
                                                                label = "",
                                                                min = 1,
                                                                max = 100,
                                                                value = 33),
                                                    helpText("Filter out the lowest 'n' percentile of the gene consensus scores from the Reference Drug Signature")
                                                    ),
                                                  
                                                    mainPanel(fluid = FALSE,
                                                              br(),
                                                              
                                                              tags$ul(
                                                                tags$li(textOutput("T2_selected_var2")),
                                                                tags$li(textOutput("T2_selected_var3"))
                                                              ),
                                                              plotlyOutput("T2_plot1"),
                                                              br(),
                                                              tableOutput("T2_value"),
                                                              br(),
                                                              br(),
                                                              br(),
                                                              downloadButton("T2_downloadData", "Download Synergy Plot Data"),
                                                              downloadButton("T2_downloadData_Drug_Sig", "Download Drug Signatures"),
                                                              DT::dataTableOutput("T2_view1")
                                                              )
                            
                                      ))
                            
                            
                            )),
                   
                  
                  
                   ### UI About ####
                   source("tabs/about.R",  local = TRUE)$value
              )
 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### 3. SERVER ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
server <- function(input, output)
                {
### Tool1  ####
source("tabs/Tool1_synergy_plot_with_reference_SERVER.R",  local = TRUE)$value
  
### Tool2  ####
  
## -> only for evaluation, delete in final version
T2_values <- reactiveValues()
                  
observeEvent(eventExpr=input$T2_Disease,ignoreInit = TRUE, {
              datasetInput_Dis <- switch(input$T2_Disease,
                                         "TCGA Bladder Urothelial Carcinoma (BLCA)" = "data/TCGA_BLCA_DE.txt", 
                                         "TCGA Breast Invasive Carcinoma (BRCA)" = "data/TCGA_BRCA_DE.txt", 
                                         "TCGA Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESC)" = "data/TCGA_CESC_DE.txt", 
                                         "TCGA Cholangiocarcinoma (CHOL)" = "data/TCGA_CHOL_DE.txt", 
                                         "TCGA Colon Adenocarcinoma (COAD)" = "data/TCGA_COAD_DE.txt", 
                                         "TCGA Esophageal Carcinoma (ESCA)" = "data/TCGA_ESCA_DE.txt", 
                                         "TCGA Glioblastoma Multiforme (GBM)" = "data/TCGA_GBM_DE.txt", 
                                         "TCGA Head and Neck Squamous Cell Carcinoma (HNSC)" = "data/TCGA_HNSC_DE.txt", 
                                         "TCGA Kidney Chromophobe (KICH)" = "data/TCGA_KICH_DE.txt", 
                                         "TCGA Kidney Renal Clear Cell Carcinoma (KIRC)" = "data/TCGA_KIRC_DE.txt", 
                                         "TCGA Kidney Renal Papillary Cell Carcinoma (KIRP)" = "data/TCGA_KIRP_DE.txt", 
                                         "TCGA Liver Hepatocellular Carcinoma (LIHC)" = "data/TCGA_LIHC_DE.txt", 
                                         "TCGA Lung Adenocarcinoma (LUAD)" = "data/TCGA_LUAD_DE.txt", 
                                         "TCGA Lung Squamous Cell Carcinoma (LUSC)" = "data/TCGA_LUSC_DE.txt", 
                                         "TCGA Pheochromocytoma and Paraganglioma (PCPG)" = "data/TCGA_PCPG_DE.txt", 
                                         "TCGA Prostate Adenocarcinoma (PRAD)" = "data/TCGA_PRAD_DE.txt", 
                                         "TCGA Rectum Adenocarcinoma (READ)" = "data/TCGA_READ_DE.txt", 
                                         "TCGA Stomach Adenocarcinoma (STAD)" = "data/TCGA_STAD_DE.txt", 
                                         "TCGA Thyroid Carcinoma (THCA)" = "data/TCGA_THCA_DE.txt", 
                                         "TCGA Uterine Corpus Endometrial Carcinoma (UCEC)" = "data/TCGA_UCEC_DE.txt", 
                                         "TCGA Glioblastoma signature (GBM.sig)" = "data/TCGA_GBM_Signature.txt", 
                                         "TCGA Colon signature (CRC.sig)" = "data/TCGA_CRC_Signature.txt", 
                                         "TCGA Breast signature (BRCA.sig)" = "data/TCGA_BRCA_Signature.txt", 
                                         "PDX GBM Group 1 (PDX.GBM1)" = "data/Table_G1_1_PDX_Group1_L1000_only.txt", 
                                         "PDX GBM Group 2 (PDX.GBM2)" = "data/Table_G2_1_PDX_Group2_L1000_only.txt", 
                                         "PDX GBM Group 3 (PDX.GBM3)" = "data/Table_G3_1_PDX_Group3_L1000_only.txt", 
                                         "PDX GBM Group 4 (PDX.GBM4)" = "data/Table_G4_1_PDX_Group4_L1000_only.txt")

              T2_Temp1 <- read.csv(file=datasetInput_Dis, header = TRUE, sep = "")
              T2_values$Disease_Signature <- T2_Temp1
              T2_values$Disease_input <- input$T2_Disease
})
                  
                  
observeEvent(eventExpr=input$T2_go_disease_signature,ignoreInit = TRUE, {
  T2_Temp2 <- read.table(text=input$T2_Disease_Text_Box,sep=",")
  colnames(T2_Temp2) <- c("Genes","log2FoldChange")
  T2_values$Disease_Signature <- T2_Temp2
  T2_values$Disease_input <- "Custom Disease Signature"
})




T2_Drug_Signatures <- reactive({
  Drugs_Sigs <- read.table(file="data/matPH3_2_1_0.2_0.3_L1000_Batch2017_Regina_removed.txt",sep ="\t", header=TRUE)
  Drugs_Sigs <- na.omit(Drugs_Sigs)
  row.names(Drugs_Sigs) <-  as.character(Drugs_Sigs$Genes)
  Drugs_Sigs <- Drugs_Sigs[,-1]
  Drugs_Sigs
                    
                  })
        
### TEXT Output: Disease Signature Gene Number ### ### ### ### ###
output$T2_selected_var2 <- renderText({
  cols <- as.character(colnames(T2_Drug_Signatures()))
  TCGA_Sig2 <- T2_values$Disease_Signature
  length2 <- intersect(cols,as.character(TCGA_Sig2$Genes))
  paste(T2_values$Disease_input,"with" ,length(length2)," genes that are measured by the L1000 Platform")
})
                  
                  
  
  
  
                }






shinyApp(ui, server)