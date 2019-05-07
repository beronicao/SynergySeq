### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### 3.1 Single Cell Sequencing #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


### Loading Single Cell Files
T3_SC_Data <- reactive({
  if(is.null(input$T3_Single_Cell_File2))
  {
    switch(input$T3_Single_Cell_File,
           "Glioblastoma" = "data/GBM_Single_Cell_only_L1000_Genes_with_coord_rounded.txt") 
  }
  else
  {
    input$T3_Single_Cell_File2$datapath
  }
})


T3_SC_Data2 <- reactive({ read.table(file=T3_SC_Data(),header = TRUE,sep = "\t",row.names=1)})


#T3_SC_Data2 <- read.table(file="data/GBM_Single_Cell_only_L1000_Genes_with_coord_rounded.txt",header = TRUE,sep = "\t",row.names=1)

### Download Template ##
output$T3_Single_Cell_File_Example <- downloadHandler(
  filename = function() {
    paste("Single_Cell_RNSAseq_File_Example", ".txt", sep = "")
  },
  content = function(file) {
    write.table(T3_SC_Data2() , file, row.names = FALSE,sep="\t")
  },
  contentType="text/plain"
) 

T3_Flavors <- read.table(file="data/Final_TCS_Signature_Correlation_02_14_2019.txt",sep="\t")
T3_Flavor_Metadata <- read.table(file="data/Final_TCS_Signature_Correlation_Metadata_02_14_2019.txt",sep="\t")
T3_Drugs <- unique(as.character(T3_Flavor_Metadata$Drug))

### Choose a Compound
output$T3_ui <- renderUI({ selectInput("T3_compound_selected", "Choose a Perturbation to Overlap",choices = T3_Drugs,selected = "JQ1-(+)")})



output$T3_plot1 <- renderPlotly({
  
  T3_IDs <- as.character(T3_Flavor_Metadata[which(T3_Flavor_Metadata$Drug==input$T3_compound_selected),"FL_IDs"])
  ##T3_IDs <- as.character(T3_Flavor_Metadata[which(T3_Flavor_Metadata$Drug=="JQ1-(+)"),"FL_IDs"])
  
  T3_Expression <- T3_SC_Data2()[,3:dim(T3_SC_Data2())[2]]
  #T3_Expression <- T3_SC_Data2[,3:dim(T3_SC_Data2)[2]]
  
  T3_Coordinates <- T3_SC_Data2()[,1:2]
  
  T3_Sigs <- T3_Flavors[T3_IDs,]
  T3_Common_Genes <- intersect(colnames(T3_Sigs),colnames(T3_Expression))
  T3_Sigs <- T3_Sigs[,T3_Common_Genes]
  T3_Expression <- T3_Expression[,T3_Common_Genes]
  
  T3_Cor1 <- t(cor(t(T3_Sigs),t(T3_Expression),method="spearman"))
  dim(T3_Cor1)
  H <- T3_Cor1[1:10,]
  
  #find the most anticorrelated flavor
  T3_Cor_metric <- apply(T3_Cor1,2,mean)
  T3_Cor_metric2 <- names(T3_Cor_metric[order(T3_Cor_metric)])[1]
  T3_Cor2 <- T3_Cor1[,T3_Cor_metric2]
  T3_Coordinates2 <- merge(T3_Coordinates,T3_Cor1[,T3_Cor_metric2],by="row.names")
  
  plot_ly(T3_Coordinates2, x = ~X_Axis, y = ~Y_Axis,color = ~y,type = 'scatter',alpha=0.7,marker = list(size = 6),colors="RdYlBu",
          mode = 'markers') %>% layout(dragmode = "select")
  
})





output$T3_plot2 <- renderPlotly({

  
  
  
  })
















