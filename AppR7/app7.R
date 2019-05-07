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
library(dplyr)
SM_MOA <- read.csv(file="data/L1000_SM_MOA2.csv", header=TRUE)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### 2. UI ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
ui <- navbarPage(  theme = shinytheme("flatly"),responsive = TRUE,
                   title = div(
                              img(src = "Logos_3.png",height="30",align="left",position="fixed",hspace="10"),
                              "SynergySeq"
                              ),
                  
                   ### UI Introduction ####
                   source("tabs/introduction_UI.R",  local = TRUE)$value,
                  
                   
                   tabPanel("Tools", tabsetPanel(
                     
                            ### UI Tool1 ####
                            source("tabs/Tool1_synergy_plot_with_reference_UI.R",  local = TRUE)$value,
                            
                            
                            ### UI Tool2 ####
                            source("tabs/Tool2_synergy_plot_no_reference_UI.R",  local = TRUE)$value,
                            
                            
                            
                            ### UI Tool3 ####
                            source("tabs/Tool3_Single_Cell_UI.R",  local = TRUE)$value
                            
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
                                         "Glioblastoma TCGA (GBM)" = "data/TCGA_GBM_Signature.txt",
                                         "Colon TCGA (CRC)" = "data/TCGA_CRC_Signature.txt",
                                         "Breast TCGA (BRCA)" = "data/TCGA_BRCA_Signature.txt",
                                         "PDX GBM Group 1" = "data/Table_G1_1_PDX_Group1_L1000_only.txt",
                                         "PDX GBM Group 2" = "data/Table_G2_1_PDX_Group2_L1000_only.txt",
                                         "PDX GBM Group 3" = "data/Table_G3_1_PDX_Group3_L1000_only.txt" ,
                                         "PDX GBM Group 4" = "data/Table_G4_1_PDX_Group4_L1000_only.txt")
              T2_Temp1 <- read.table(file=datasetInput_Dis,header = TRUE,sep = "\t")
              T2_values$Disease_Signature <- T2_Temp1
              T2_values$Disease_input <- input$T2_Disease
})
                  
                  
observeEvent(eventExpr=input$T2_go_disease_signature,ignoreInit = TRUE, {
  T2_Temp2 <- read.table(text=input$T2_Disease_Text_Box,sep=",")
  colnames(T2_Temp2) <- c("Genes","log2FoldChange")
  T2_values$Disease_Signature <- T2_Temp2
  T2_values$Disease_input <- "Custom Disease Signature"
})


T2_TCS_Flavors <- read.delim(file="data/Final_TCS_Signature_Correlation_02_14_2019.txt")
T2_TCS_Flavor_Metadata <- read.delim(file="data/Final_TCS_Signature_Correlation_Metadata_02_14_2019.txt")
T2_Similarity <- read.table(file="data/TCS_Weighted_SR.txt",row.names=1,header=TRUE,check.names = FALSE)
T2_Similarity <- as.matrix(T2_Similarity)

T2_Table_Final <- reactive({
  
T2_Disease_Sig <- T2_values$Disease_Signature
T2_Disease_Sig <- read.delim(file="data/TCGA_GBM_Signature.txt")
T2_Common_Genes <- intersect(as.character(T2_Disease_Sig$Genes),colnames(T2_TCS_Flavors))
T2_Disease_Sig2 <- T2_Disease_Sig[which(T2_Disease_Sig$Genes %in% T2_Common_Genes),]
T2_Disease_Sig3 <- as.numeric(as.character(T2_Disease_Sig2$log2FoldChange)) 
names(T2_Disease_Sig3) <- as.character(T2_Disease_Sig2$Genes)
T2_TCS_Flavors2 <- T2_TCS_Flavors[,names(T2_Disease_Sig3)]
T2_TCS_Flavors2[abs(T2_TCS_Flavors2)==1] <- 0
T2_Drugs <- unique(as.character(T2_TCS_Flavor_Metadata$Drug))


### STEP1. Calculate Individual Drug Discordance Ratio 
T2_Disc_Ratio <- as.data.frame(apply(T2_TCS_Flavors2,1,function(t2_x){
  
  t2_prod <- t2_x*T2_Disease_Sig3
  t2_a <- abs(length(t2_prod[t2_prod<0]))
  t2_b <- abs(length(t2_prod[t2_prod>0]))
  if(t2_b==0){t2_b<-1}
  t2_c <- t2_a/t2_b
}))
colnames(T2_Disc_Ratio) <- "DR"
T2_Disc_Ratio2 <- merge(T2_Disc_Ratio,T2_TCS_Flavor_Metadata,by.x="row.names",by.y="FL_IDs")


### STEP2. Filter for Drugs that have a DR > 2
T2_Drugs2 <- unique(as.character(T2_Disc_Ratio2[which(T2_Disc_Ratio2$DR>=2),"Drug"]))

I <- intersect(colnames(T2_Similarity),row.names(T2_Similarity))

T2_Drugs2 <- intersect(T2_Drugs2,I)

# T3 <- T2_Disc_Ratio2 %>% group_by(Drug) %>% summarise(
#   min = min(DR),
#   max = max(DR),
#   mean =mean(DR))
# ggplot(T2_Disc_Ratio2, aes(x = Drug, y = DR)) +
#   geom_boxplot()

### STEP3. Calculate the max of the DR Flavors
T2_Table1 <- data.frame()
for (t2_drug in T2_Drugs2)
{
  t2_1 <- T2_Disc_Ratio2[which(T2_Disc_Ratio2$Drug==t2_drug),]
  t2_2 <- t2_1[order(t2_1$DR,decreasing=TRUE),]
  t2_3 <- as.character(t2_2[1,"Row.names"])
  t2_4 <- as.numeric(as.character(t2_2[1,"DR"]))
  t2_vec <- c(t2_drug,t2_3,t2_4)
  T2_Table1 <- rbind(T2_Table1,t2_vec,stringsAsFactors=FALSE)
}

T2_Table1 <- T2_Disc_Ratio2[which(T2_Disc_Ratio2$DR >=2),c("Drug","Row.names","DR")]
colnames(T2_Table1) <- c("Drug","Flavor","DR")


### NOT USED STEP3. Use the Highest DR among the flavors
# for (t2_drug in T2_Drugs2)
# {
#   t2_1 <- T2_Disc_Ratio2[which(T2_Disc_Ratio2$Drug==t2_drug),]
#   t2_2 <- t2_1[order(t2_1$DR,decreasing=TRUE),]
#   t2_3 <- as.character(t2_2[1,"Row.names"])
#   t2_4 <- as.numeric(as.character(t2_2[1,"DR"]))
#   t2_vec <- c(t2_drug,t2_3,t2_4)
#   T2_Table1 <- rbind(T2_Table1,t2_vec,stringsAsFactors=FALSE)
# }
# colnames(T2_Table1) <- c("Drug","Flavor","DR")




### Check how many unique
t2_comb1 <- t(combn(as.character(T2_Table1$Flavor),2))
T2_TCS_Flavors3 <- as.matrix(T2_TCS_Flavors2)
T2_TCS_Flavors3[abs(T2_TCS_Flavors3)==1] <-0
#t2_fl_1 <- "alisertib$FL1"
#t2_fl_2 <- "GSK-1070916$FL2"

#Start.time <- Sys.time()
# T2_Table2 <- apply(t2_comb1,1,function(t2_x){
# 
#   t2_fl_1 <- as.character(unlist(t2_x[1]))
#   t2_fl_2 <- as.character(unlist(t2_x[2]))
#   t2_dr_1 <- as.character(T2_TCS_Flavor_Metadata[which(T2_TCS_Flavor_Metadata$FL_IDs==t2_fl_1),"Drug"])
#   t2_dr_2 <- as.character(T2_TCS_Flavor_Metadata[which(T2_TCS_Flavor_Metadata$FL_IDs==t2_fl_2),"Drug"])
#   if (t2_dr_1 %in% row.names(T2_Similarity) & t2_dr_2 %in% colnames(T2_Similarity))
#   {
#   t2_sim <- T2_Similarity[t2_dr_1,t2_dr_2]
#   t2_sig1 <- unlist(T2_TCS_Flavors3[t2_fl_1,])
#   t2_sig2 <- unlist(T2_TCS_Flavors3[t2_fl_2,])
#   t2_p <- t2_sig1*t2_sig2
#   
#   #Drug1
#   t2_sig1[names(t2_p[t2_p<0])] <- 0
#   ts_drug1 <- T2_Disease_Sig3*t2_sig1
#   t2_drug1_a <- length(ts_drug1[ts_drug1<0])
#   #t2_drug1_a <- sum(abs(T2_Disease_Sig3[ts_drug1<0]))
#   t2_drug1_b <- length(ts_drug1[ts_drug1>0])
#   #t2_drug1_b <- sum(abs(T2_Disease_Sig3[ts_drug1>0]))
#   if (t2_drug1_b==0) { t2_drug1_b<-1}
#   t2_drug1_c <- t2_drug1_a/t2_drug1_b
#   t2_drug1_d <- t2_drug1_c *t2_drug1_a
#   
#   #Drug2
#   t2_sig2[names(t2_p[t2_p<0])] <- 0
#   ts_drug2 <- T2_Disease_Sig3*t2_sig2
#   t2_drug2_a <- length(ts_drug2[ts_drug2<0])
#   #t2_drug2_a <- sum(abs(T2_Disease_Sig3[ts_drug2<0]))
#   t2_drug2_b <- length(ts_drug2[ts_drug2>0])
#   #t2_drug2_b <- sum(abs(T2_Disease_Sig3[ts_drug2>0]))
#   if (t2_drug2_b==0) { t2_drug2_b<-1}
#   t2_drug2_c <- t2_drug2_a/t2_drug2_b
#   t2_drug2_d <- t2_drug2_c * t2_drug2_a
#   t2_sum <- t2_drug2_a+t2_drug2_a
#   
#   #Combined Drug
#   t2_sig3 <- t2_sig1 +t2_sig2
#   t2_drug3_a <- length(t2_sig3[t2_sig3<0])
#   
#   t2_mean <- mean(c(t2_drug1_c,t2_drug2_c))
#   t2_mean2 <- mean(c(t2_drug1_d,t2_drug2_d))
#   if(t2_drug1_c>2 & t2_drug2_c > 2 & t2_drug1_a>20 & t2_drug2_a>20 )
#   {
#     return(c(t2_fl_1,t2_fl_2,t2_sim,t2_drug1_c,t2_drug2_c,t2_mean,t2_drug1_a,t2_drug2_a,t2_sum,t2_drug1_d,t2_drug2_d,t2_mean2))
# 
#   }
#   }
# })
# #End.time <- Sys.time()
# #End.time-Start.time
# T2_Table3 <-Filter(length, T2_Table2)
# T2_Table4 <- data.frame(t(sapply(T2_Table3,c)))
# T2_Table4$X1 <- as.character(T2_Table4$X1)
# T2_Table4$X2 <- as.character(T2_Table4$X2)
# T2_Table4$X3 <- as.numeric(as.character(T2_Table4$X3))
# T2_Table4$X4 <- as.numeric(as.character(T2_Table4$X4))
# T2_Table4$X8 <- as.numeric(as.character(T2_Table4$X8))
# T2_Table4$X5 <- as.numeric(as.character(T2_Table4$X5))
# T2_Table4$X6 <- as.numeric(as.character(T2_Table4$X6))
# T2_Table4$X9 <- as.numeric(as.character(T2_Table4$X9))
# T2_Table4$X10 <- as.numeric(as.character(T2_Table4$X10))
# T2_Table4$X11 <- as.numeric(as.character(T2_Table4$X11))
# T2_Table4$X12 <- as.numeric(as.character(T2_Table4$X12))

#write.table(T2_Table4,file="data/T2_Table4.txt",sep="\t")

#plot_ly(T2_Table4, x = ~X3, y = ~X12,type = 'scatter',alpha=0.7,marker = list(size = 10),
#        mode = 'markers',hoverinfo= 'text',text=~paste(X1," ",X2)) %>% layout(dragmode = "select")

# T2_FL_ID1 <- unique(c(as.character(T2_Table4$X1),as.character(T2_Table4$X2)))
# T2_Table4$ID <- paste(as.character(T2_Table4$X1),as.character(T2_Table4$X2),sep="")
# T2_Table4a <- T2_Table4[,c("ID","X1","X4")]
# T2_Table4a$Col <- "X1"
# colnames(T2_Table4a) <- c("ID","Drug","DR","Col")
# T2_Table4b <- T2_Table4[,c("ID","X2","X5")]
# T2_Table4b$Col <- "X2"
# colnames(T2_Table4b) <- c("ID","Drug","DR","Col")
# T2_Table5 <- rbind(T2_Table4a,T2_Table4b)
# 
# T2_Final2 <- data.frame()
# for (t2_id in T2_FL_ID1)
# {
# 
#   t2_vector1 <-T2_Table5[which(T2_Table5$Drug==t2_id),]
#   t2_vector2 <- as.numeric(as.character(t2_vector1$DR))
#   t2_vector3 <- scale(t2_vector2)
# 
#   #using POMS instead of z score https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4569815/#B14
#   t2_vector1$Z <- t2_vector3
#   t2_vector1$P <-(t2_vector2-min(t2_vector2))/(max(t2_vector2)-min(t2_vector2))
#   t2_vector1$N <- t2_vector2*t2_vector1$P
#   T2_Final2 <- rbind(T2_Final2,t2_vector1,stringsAsFactors=FALSE)
# 
# 
# }
# library(reshape)
# T2_Final3 <- dcast(T2_Final2,ID~Col,value.var="N",fun.aggregate=mean)
# 
# T2_Final3 <- cast(T2_Final2,ID~Drug,value.var=N)
# T2_Final3 <- aggregate(T2_Final2$P,by=list(as.character(T2_Final2$ID)),FUN=mean)
# T2_Final4 <- T2_Final2[,c("ID","Col","Z")]
# 
# 
# 
# 
# tq <- T2_FL_ID1[2]
# for (tq in T2_FL_ID1)
# {
#   
#   T2_Table5a <- T2_Table4[which(T2_Table4$X1 ==tq | T2_Table4$X2 ==tq),c("X1","X2","X6")]
#   
# }
# 
# T2_Table4





})



output$T2_plot1 <- renderPlotly({
  
  T2_Table_Final <- read.table(file="data/T2_Table4.txt")
  plot_ly(T2_Table_Final, x = ~X3, y = ~X6,type = 'scatter',alpha=0.7,marker = list(size = 10),
          mode = 'markers',hoverinfo= 'text',text=~paste(X1," ",X2,'<br>',X3," ",X4)) %>% layout(dragmode = "select")
})



#NOT USED
# 
# ### IMPORTANT 04/19/2019### <<<<<<<<<<<<<<<<<<<<<<<<<<<<
# t2_id <- "ACY-1215$FL1"
# T2_FL_ID1 <- unique(as.character(T2_Table5$Drug))
# 
# 
# 
# 
# 
# T2_Final5 <- dcast(T2_Final4,ID~Col,fun.aggregate=mean)
# 
# 
# T2_Final6 <- apply(T2_Final5,1, function(x){
#  m <- c(as.numeric(as.character(unlist(x[2]))),as.numeric(as.character(unlist(x[3]))))
#   mean(m)
#   
#   
# })
# T2_Final5$Mean <- T2_Final6
# 
# 
# 
# 
# 
# 
# 
# 
# 
#   t2_union_table[,names(t2_p[t2_p<0])] <-0
#   
#   t2_gene1 <- names(t2_sig1[t2_sig1!=0])
#   t2_gene2 <- names(t2_sig2[t2_sig2!=0])
#   t2_diff1 <- length(setdiff(t2_gene1,t2_gene2))
#   t2_diff2 <- length(setdiff(t2_gene2,t2_gene1))
#   
#   t2_union1 <- union(t2_gene1,t2_gene2)
#   
#   
#   
#   t2_p <- t2_sig1[t2_union1]*t2_sig2[t2_union1]
#   names(t2_p[t2_p>0])
#   t2_union_table[,names(t2_p[t2_p<0])] <-0
#   ts_drug1 <- T2_Disease_Sig3*
#  
#   
#   # Remove genes that are going opposite directions
#   t2_x <- t2_union_table[,1]
#   t2_union_table2 <- apply(t2_union_table,2,function(t2_x){
#     
#     t2_p <- t2_x[1]*t2_x[2]
#     names(t2_p[t2_p<0])
#     if(t2_p < 0){
#       t2_final <-0
#     }
#     
#     
#     if(t2_p >= 0){
#       
#       if(t2_x[1] > 0)
#       {
#         t2_final <- max( t2_x[1],t2_x[2])
#         
#       }
#       
#       
#       if(t2_x[1] <0)
#       {
#         t2_final <- min( t2_x[1],t2_x[2])
#         
#       }
#       
#       if(t2_x[1] == 0)
#       {
#         t2_final <- t2_x[2]
#         
#       }
#     }
#     
#     t2_final
#     
#   })
#   
#   names(t2_union_table2) <- t2_union1
#   
#   T2_Disease_Sig4 <- T2_Disease_Sig3[t2_union1]
#   t2_p2 <- t2_union_table2*T2_Disease_Sig4
#   
#   t2_a <- length(t2_p2[which(t2_p2<0)])
#   t2_b <- length(t2_p2[which(t2_p2>0)])
#   if(t2_b==0){t2_b<-1}
#   t2_c <- t2_a/t2_b
#   t2_vect <- c(t2_fl_1,t2_fl_2,t2_c)
#   
# 
# End.time <- Sys.time()
# 
# T2_Table3 <- t(T2_Table2)
# 
# T2_Table3 <-Filter(length, T2_Table2)
# 
# T2_Table4 <- data.frame(t(sapply(T2_Table3,c)))
# 












# T2_Drug_Sig_Threshold <- 50



# T2_Drug_Pairs <- reactive({
#   
#   ##T2_Drugs <- row.names(Drugs_Sigs)
#   T2_Drugs <- row.names(T2_Drug_Signatures())
#   T2_Drug_Combos <- t(combn(T2_Drugs,2))
# })  
#   

# T2_Table1 <- reactive({
#   T2_Drug_SR <- read.table(file=T2_Drug_Similarity_Metrics(),sep="\t",header=TRUE,check.names = FALSE)
#   ##T2_Drug_SR <- read.table(file="data/TCS_Weighted_SR.txt",sep="\t",header=TRUE,check.names = FALSE)
#   T2_rows <- as.character(T2_Drug_SR$DrugA)
#   row.names(T2_Drug_SR) <- T2_rows
#   T2_Drug_SR <- T2_Drug_SR[,-1]
#   T2_Drug_SR <- as.matrix(T2_Drug_SR)
#   H <- T2_Drug_SR[1:10,1:10]
#   T2_Drug_SR2 <- melt(T2_Drug_SR)
#   
#   #T2_Drug_DR <- 
#   
#   
# })


### Text Output: Disease Signature Gene Number ###
output$T2_selected_var2 <- renderText({
  cols <- as.character(colnames(T2_TCS_Flavors))
  TCGA_Sig2 <- T2_values$Disease_Signature
  length2 <- intersect(cols,as.character(TCGA_Sig2$Genes))
  paste(T2_values$Disease_input,"with" ,length(length2)," genes that are measured by the L1000 Platform")
})




### Tool3  ####
source("tabs/Tool3_Single_Cell_SERVER.R",  local = TRUE)$value


                  

                }






shinyApp(ui, server)
