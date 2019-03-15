### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### 2.1 Synergy Plot #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

values <- reactiveValues()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.1 Loading Disease Signatures  
datasetInput_Dis <- reactive({
  switch(input$disease,
         "Glioblastoma TCGA (GBM)" = "data/TCGA_GBM_Signature.txt",
         "Colon TCGA (CRC)" = "data/TCGA_CRC_Signature.txt",
         "Breast TCGA (BRCA)" = "data/TCGA_BRCA_Signature.txt",
         "PDX GBM Group 1" = "data/Table_G1_1_PDX_Group1_L1000_only.txt",
         "PDX GBM Group 2" = "data/Table_G2_1_PDX_Group2_L1000_only.txt",
         "PDX GBM Group 3" = "data/Table_G3_1_PDX_Group3_L1000_only.txt" ,
         "PDX GBM Group 4" = "data/Table_G4_1_PDX_Group4_L1000_only.txt") 
})  

TCGA_Sig <- reactive({ read.table(file=datasetInput_Dis(),header = TRUE,sep = "\t")})



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
### 2.1.6 Calculating Similarity and Discordance Scores
Final3 <- reactive({
  JQ1_Sig_v <- JQ1_Sig()
  #JQ1_Sig_v <- t(Drugs_Sigs["JQ1",])
  JQ1_Sig_v<- data.frame(JQ1_Sig_v,Genes=row.names(JQ1_Sig_v))
  colnames(JQ1_Sig_v) <- c("JQ1","Genes")
  TCGA_Sig_v <- TCGA_Sig()
  #TCGA_Sig_v <- read.table(file="data/TCGA_GBM_Signature.txt",header = TRUE,sep = "\t")
  V <- unique(as.numeric(as.character(JQ1_Sig_v$JQ1)))
  V2 <- max(abs(V))
  V4 <- as.numeric(input$bins)/100
  #V4 <- as.numeric(30/100)
  V5 <- V2*V4
  V6 <- round(V5)
  JQ1_Sig_v$JQ1 <- as.numeric(as.character(JQ1_Sig_v$JQ1))
  JQ1_Sig_v2 <- JQ1_Sig_v[which(JQ1_Sig_v$JQ1 > V6 | JQ1_Sig_v$JQ1 < -V6),]
  JQ1_Sig_v2 <- JQ1_Sig_v2[which(JQ1_Sig_v2$JQ1 > 0 | JQ1_Sig_v2$JQ1 < 0),]
  jq1_genes <- as.character(JQ1_Sig_v2$Genes)
  Drugs_Sigs2 <- Drugs_SigsR()[,jq1_genes]
  tt <-  apply(Drugs_Sigs2,1,function(x){x*JQ1_Sig_v2$JQ1})
  
  #this is to calculate how similar are the drug signatures to jq1 (ratio of # of genes that have same direction with the JQ1 signature devided with # of genes that are disocrdant to JQ1)
  tt2 <-  apply(tt,2,function(x) { a <- sum(x>0)
  b <- sum(x<0)
  b[b==0] <- 1 # replace b with 1 so we can devide by that number. Another way would be to add one to both a and b
  c <- a/b
  c
  })
  tt3 <- as.data.frame(tt2)
  tmax <- max(tt3)
  tt4 <- tt3/tmax
  
  #this is to calculate the ratio of the genes that are discordant to the Disease Signature (and are not affected by JQ1) devided by the # of genes that are concordant to the Disease Signature (and non-JQ1)
  jq1_genes2 <- as.character(JQ1_Sig_v[which(JQ1_Sig_v$JQ1==0),2])
  #01/23/2019 jq1_genes2 <- as.character(JQ1_Sig_v[which(JQ1_Sig_v$JQ1 < V6 & JQ1_Sig_v$JQ1 > -V6),2])
  
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
    aa[aa==0] <- 1 # replace bb with 1 so we can devide by that number. Another way would be to add one to both aa and bb
    cc <- bb/aa
    
    dd<- c(aa,bb,cc)
  })
  pp3 <- as.data.frame(t(pp2))
  pmax <- max(pp3[,3])
  pp4 <- pp3/pmax
  Final <- merge(pp4,tt4,by="row.names")
  #text2 <- paste(input$signature,"_Orthogonality",sep="")
  colnames(Final) <- c("Drug","Disease_Same","Disease_Opp","Disease_Discordance","Reference_Drug_Orthogonality")
  Final
  #values$Final2 <- data.frame(Final)
  Final[,2] <- round(Final[,2],digits=3)
  Final[,3] <- round(Final[,3],digits=3)
  
  values$Final2 <- merge(Final,SM_MOA,by.x="Drug",by.y="Drugs",all.x = TRUE)
  
})


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2.1.7 Synergy Plot
output$plot1 <- renderPlotly({
  plot_ly(Final3(), x = ~Reference_Drug_Orthogonality, y = ~Disease_Discordance,type = 'scatter',alpha=0.7,marker = list(size = 14),
          mode = 'markers',hoverinfo= 'text',text=~paste(Drug,'<br>',"Ratio:",Disease_Discordance,Reference_Drug_Orthogonality)) %>% layout(dragmode = "select")
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
