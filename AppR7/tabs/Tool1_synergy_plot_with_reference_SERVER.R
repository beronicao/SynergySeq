### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### 2.1 Synergy Plot #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

values <- reactiveValues()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.1 Loading Disease Signatures  
datasetInput_Dis <- reactive({
  switch(input$disease,
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
})  

TCGA_Sig <- reactive({ read.csv(file=datasetInput_Dis(),sep = "")})

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.1 Loading Drug Signatures

datasetInput_Dataset <- reactive({
  if (is.null(input$L1000_Dataset))
  {return()}
  switch(input$L1000_Dataset,
         "LINCS L1000 Dec 2015" = "data/OUT3_noDMSO_noUntreated_Regina_removed.txt",
         "LINCS L1000 March 2017" = "data/matPH3_2_1_0.2_0.3_L1000_Batch2017_Regina_removed.txt"
  )
})

Drugs_SigsR <- reactive({
  Drugs_Sigs <- read.table(file=datasetInput_Dataset(),sep ="\t", header=TRUE)
  Drugs_Sigs <- na.omit(Drugs_Sigs)
  row.names(Drugs_Sigs) <-  as.character(Drugs_Sigs$Genes)
  Drugs_Sigs <- Drugs_Sigs[,-1]
  Drugs_Sigs
})

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.3 List of Drug Names

Drugs <- reactive({ sort(row.names(Drugs_SigsR()))})

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.4 Selecting a Reference Drug

output$ui <- renderUI({ selectInput("signature", "Select a Reference Drug",choices = Drugs(),selected = "GBM_JQ1")})

JQ1_Sig <- reactive({
  T <-t(Drugs_SigsR()[as.character(input$signature),])
  
})

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.5 Reference Drug Threshold

output$selected_var <- renderText({
  value2 <- as.numeric(100 - as.numeric(input$bins))
  paste("You have selected to use a ", value2,"% threshold for the gene consensus score. Genes with the lowest ",input$bins,"% scores will be filtered out",sep = "")
})

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.6 Selecting a Measure of Bioactivity

datasetInput_Bioactivities <- reactive({
  if (is.null(input$bioactivities))
  {return()}
  
  switch(input$bioactivities,
         "IC50" = "data/Bioactivities_IC50.txt",
         "Kd" = "data/Bioactivities_Kd.txt",
         "Ki" = "data/Bioactivities_Ki.txt",
         "Potency" = "data/Bioactivities_Potency.txt"
  )
})

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.6 Selecting a Measure of Bioactivity

Drugs_TargsR <- reactive({
  Drugs_Targs <- read.table(file=datasetInput_Bioactivities(), sep ="\t", header=TRUE)
  Drugs_Targs.moa <- merge(Drugs_Targs, SM_MOA, by = "Drugs")
  Drugs_Targs.moa
})

Targets <- reactive({ sort(as.character(Drugs_TargsR()[,'target_gene_symbol']))})

output$ui_bioact <- renderUI({ selectInput("target", "Select a molecular target", choices = Targets())})

targetChoice <- reactive({
  as.character(input$target)
})

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.8 Calculating Similarity and Discordance Scores

Final3 <- reactive({
  JQ1_Sig_v <- JQ1_Sig()
  JQ1_Sig_v<- data.frame(JQ1_Sig_v,Genes=row.names(JQ1_Sig_v))
  colnames(JQ1_Sig_v) <- c("JQ1","Genes")
  TCGA_Sig_v <- TCGA_Sig()
  V <- unique(as.numeric(as.character(JQ1_Sig_v$JQ1)))
  V2 <- max(abs(V))
  V4 <- as.numeric(input$bins)/100
  V5 <- V2*V4
  V6 <- round(V5)
  JQ1_Sig_v$JQ1 <- as.numeric(as.character(JQ1_Sig_v$JQ1))
  JQ1_Sig_v2 <- JQ1_Sig_v[which(JQ1_Sig_v$JQ1 > V6 | JQ1_Sig_v$JQ1 < -V6),]
  JQ1_Sig_v2 <- JQ1_Sig_v2[which(JQ1_Sig_v2$JQ1 > 0 | JQ1_Sig_v2$JQ1 < 0),]
  jq1_genes <- as.character(JQ1_Sig_v2$Genes)
  Drugs_Sigs2 <- Drugs_SigsR()[,jq1_genes]
  tt <-  apply(Drugs_Sigs2,1,function(x){x*JQ1_Sig_v2$JQ1})
  
  # Calculate how similar the drug signatures are to jq1 (ratio of # of genes that have same direction with the JQ1 signature divided by # of genes that are disocrdant to JQ1)
  tt2 <-  apply(tt,2,function(x) { 
    a <- sum(x>0)
    b <- sum(x<0)
    b[b==0] <- 1 # replace b with 1 so we can devide by that number. Another way would be to add one to both a and b
    c <- a/b
    # c
    return(c(a,b,c))
  })
  
  tt3 <- as.data.frame(t(tt2))
  
  # tt3 <- as.data.frame(tt2)
  # tmax <- max(tt3)
  # tt4 <- tt3/tmax
  ### ### ### ### ### ### ### ### ### ### ### ###
  
  # Calculate the ratio of the genes that are discordant to the Disease Signature (and are not affected by JQ1) divided by the # of genes that are concordant to the Disease Signature (and non-JQ1)
  jq1_genes2 <- as.character(JQ1_Sig_v[which(JQ1_Sig_v$JQ1==0),2])

  genes <- colnames(Drugs_SigsR())
  TCGA_Sig_v$log2FoldChange <- as.numeric(TCGA_Sig_v$log2FoldChange)
  TCGA_Sig_v2 <- TCGA_Sig_v[,-1,drop=FALSE]
  TCGA_Sig_v3 <- aggregate(TCGA_Sig_v2, by = list(TCGA_Sig_v$Genes),FUN=mean) # this is in case we have duplicate gene symbols in the signature
  non_jq1_genes <- intersect(as.character(TCGA_Sig_v3$Group.1),jq1_genes2)
  Drugs_Sigs3 <- Drugs_SigsR()[,non_jq1_genes]
  row.names(TCGA_Sig_v3) <- as.character(TCGA_Sig_v3$Group.1)
  TCGA_Sig_v4 <- TCGA_Sig_v3[non_jq1_genes,]
  
  pp <-  apply(Drugs_Sigs3,1,function(x){x*TCGA_Sig_v4$log2FoldChange})
  pp2 <-  apply(pp,2,function(x) { 
    aa <- sum(x>0)
    bb <- sum(x<0)
    aa[aa==0] <- 1 # replace aa with 1 so we can divide by that number. Another way would be to add one to both aa and bb
    cc <- bb/aa
    
    dd<- c(aa,bb,cc)
  })
  
  pp3 <- as.data.frame(t(pp2))
  
  # pmax <- max(pp3[,3])
  # pp4 <- pp3/pmax
  # Final <- merge(pp4,tt4,by="row.names")
  # colnames(Final) <- c("Drug","Disease_Same","Disease_Opp","Disease_Discordance","Reference_Drug_Orthogonality")
  
  Final <- merge(pp3[,3,drop=F],tt3[,3,drop=F],by="row.names")
  colnames(Final) <- c("Drug","Disease_Discordance", "Reference_Drug_Orthogonality")
  Final
  
  Final[,2] <- round(Final[,2],digits=3)
  Final[,3] <- round(Final[,3],digits=3)
  
  bioact.moa <- Drugs_TargsR()
  
  # values$Final2 <- merge(Final,SM_MOA,by.x="Drug",by.y="Drugs",all.x = TRUE)
  
  f3 <- merge(Final, bioact.moa, by.x = "Drug", by.y = "Drugs", all.x = TRUE)
  
  f3.ord <-f3[with(f3, order(-Disease_Discordance, Reference_Drug_Orthogonality)), ]
  
  target <- targetChoice()
  drugs.sameTargets1 <- na.omit(as.data.frame(f3.ord[f3.ord$target_gene_symbol==paste0(target),]))
  colnames(drugs.sameTargets1)[1] <- c("Drug")
  
  if (nrow(drugs.sameTargets1)>2){
    pal <- brewer.pal(nrow(drugs.sameTargets1),"BrBG")
  } else{
    pal <- brewer.pal(3,"BrBG")
  }
  
  values$Final2 <- merge(Final,SM_MOA,by.x="Drug",by.y="Drugs",all.x = TRUE)
  
})

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.9 Synergy Plot
output$plot1 <- renderPlotly({
  p1 <- plot_ly(Final, x = ~Reference_Drug_Orthogonality, y = ~Disease_Discordance, type = 'scatter', alpha=0.7, marker = list(size = 14),
                mode = 'markers', hoverinfo = 'text', text = ~paste(Drug,'<br>', "Ratio:", Disease_Discordance, Reference_Drug_Orthogonality)) %>% layout(dragmode = "select")
  
  p1 %>%
    add_data(drugs.sameTargets1) %>%
    add_markers(x = ~Reference_Drug_Orthogonality, y = ~Disease_Discordance, color = ~Median, colors = pal) %>% colorbar(title = paste0("Median ", input$bioactivities))
})

output$table1 <- renderTable({
  
  s <- event_data("plotly_selected")
  values$Final2$Reference_Drug_Orthogonality <- as.character(values$Final2$Reference_Drug_Orthogonality)
  values$Final2$Disease_Discordance <- as.character(values$Final2$Disease_Discordance)
  values$Final2[which(values$Final2$Reference_Drug_Orthogonality %in% s$x & values$Final2$Disease_Discordance %in% s$y)  ,]
})

output$view1 <- DT::renderDataTable({ 
  s <- event_data("plotly_selected")
  values$Final2$Reference_Drug_Orthogonality <- as.character(values$Final2$Reference_Drug_Orthogonality)
  values$Final2$Disease_Discordance <- as.character(values$Final2$Disease_Discordance)
  values$Final2[which(values$Final2$Reference_Drug_Orthogonality %in% s$x & values$Final2$Disease_Discordance %in% s$y)  ,]
  #values$Final2$Reference_Drug_Orthogonality <- round(as.numeric(values$Final2$Reference_Drug_Orthogonality),digits=2 )
  values$Final2[,c(1,2,3,6,7)]
})


output$selected_var2 <- renderText({
  cols <- as.character(colnames(Drugs_SigsR()))
  TCGA_Sig2 <- TCGA_Sig()
  length2 <- intersect(cols,as.character(TCGA_Sig2$Genes))
  paste("The Disease Signature has ",length(length2)," genes that are measured by the L1000 Platform")
})


output$selected_var3 <- renderText({
  JQ1_Sig3 <- JQ1_Sig()
  JQ1_Sig3<-  data.frame(JQ1_Sig3,Genes=row.names(JQ1_Sig3))
  colnames(JQ1_Sig3) <- c("JQ1","Genes")
  T <- unique(as.numeric(as.character(JQ1_Sig3$JQ1)))
  T2 <- max(abs(T))
  T4 <- as.numeric(input$bins)/100
  T5 <- T2*T4
  T6 <- round(T5)
  JQ1_Sig3$JQ1 <- as.numeric(as.character(JQ1_Sig3$JQ1))
  JQ1_Sig3_2 <- JQ1_Sig3[which(JQ1_Sig3$JQ1 > T6 | JQ1_Sig3$JQ1 < -T6),]
  JQ1_Sig3_2 <- JQ1_Sig3_2[which(JQ1_Sig3_2$JQ1 > 0 | JQ1_Sig3_2$JQ1 < 0),]
  paste("The",input$signature," signature has",length(JQ1_Sig3_2$Genes),"genes that are measured by the L1000 Platform")
})


output$downloadData <- downloadHandler(
  filename = function() {
    paste("output", ".txt", sep = "")
  },
  
  content = function(file) {
    write.table(values$Final2 , file, row.names = TRUE,sep="\t")
    
  },
  contentType="text/plain"
)


output$downloadData_Drug_Sig <- downloadHandler(
  filename = function() {
    paste("Drug_Signature_Table", ".txt", sep = "")
  },
  content = function(file) {
    write.table(Drugs_SigsR() , file, row.names = TRUE,sep="\t")
    
  },
  contentType="text/plain"
) 


output$downloadExample <- downloadHandler(
  filename = function() {
    paste("Drug_Signature_Example", ".txt", sep = "")
  },
  content = function(file) {
    write.table(TCGA_Sig() , file, row.names = FALSE,sep="\t")
    
  },
  contentType="text/plain"
) 
