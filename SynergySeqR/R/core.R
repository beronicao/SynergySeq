#' Drugs_SigsR
#'
#' Reads input data files and calculates similarity score of each drug to reference, ratio of genes discordant to disease signature.
#' Returns a data frame of similarity scores.
#'
#' @param datasetInput_Dataset Specify a dataset of drug signatures. (Default: LINCS L1000 Dec 2015) For custom drug signature table, input name of a tab-delimited file with drug/molecule names in the first column and gene symbols in the first row.   
#' @param refDrug Character vector specifying a reference drug signature to use from datasetInput_Dataset. (Default: JQ1)
#' @param n_bins Numeric value to filter out the lowest 'n' percentile of the gene consensus scores from the Reference Drug Signature. (Default value: 33)
#' @param datasetInput_Dis Specify a gene expression signature dataset for a disease of interest. (Default: Glioblastoma TCGA) 
#' @param ... Optional arguments passed to read.table().
#'
#' @return a data frame 
#'
#' @family Drugs_SigsR functions
#'
#' @export


Drugs_SigsR <- function(datasetInput_Dataset=NULL, refDrug=NULL, n_bins=33, datasetInput_Dis=NULL, options_Bioactiv=NULL, molecTarget=NULL, ...){
  if (is.null(datasetInput_Dataset)==TRUE){
    inputFile <- "data/OUT3_noDMSO_noUntreated_Regina_removed.txt"
    message('No input dataset specified. Defaulting to "LINCS L1000 Dec 2015". ')
  } else {
    
    if (datasetInput_Dataset=="LINCS L1000 Dec 2015"){
      inputFile <- "data/OUT3_noDMSO_noUntreated_Regina_removed.txt"
    } else {
      
      if (datasetInput_Dataset=="LINCS L1000 March 2017"){
        inputFile <- "data/matPH3_2_1_0.2_0.3_L1000_Batch2017_Regina_removed.txt"
      } else {
        inputFile <- datasetInput_Dataset
        #next
      }
    }
  }
  Drugs_Sigs <- read.table(file=inputFile, sep ="\t", header=TRUE)
  Drugs_Sigs <- na.omit(Drugs_Sigs)
  row.names(Drugs_Sigs) <-  as.character(Drugs_Sigs[,1])
  Drugs_Sigs <- Drugs_Sigs[,-1]
  # return(Drugs_Sigs)
  
  if (is.null(refDrug)==TRUE){
    JQ1_Sig <- t(Drugs_Sigs[1,])
    message(paste0('No reference drug signature specified. Defaulting to ', colnames(JQ1_Sig)))
  } else {
    
    JQ1_Sig <- t(Drugs_Sigs[paste0(refDrug),])
    
  }
  
  value2 <- as.numeric(100 - as.numeric(n_bins))
  message(paste("Using a ", value2,"% threshold for the gene consensus score. Genes with the lowest ",n_bins,"% scores will be filtered out",sep = ""))
  
  
  if (is.null(datasetInput_Dis)==TRUE){
    input_disease <- "GBM.sig"
    message('No input disease specified. Defaulting to Glioblastoma TCGA (GBM.sig). ')
  } else {
    
    input_disease <- datasetInput_Dis
  }
  
  input_disease <- switch(input_disease,
                          "BLCA" = "data/TCGA_BLCA_DE.txt", 
                          "BRCA" = "data/TCGA_BRCA_DE.txt", 
                          "CESC" = "data/TCGA_CESC_DE.txt", 
                          "CHOL" = "data/TCGA_CHOL_DE.txt", 
                          "COAD" = "data/TCGA_COAD_DE.txt", 
                          "ESCA" = "data/TCGA_ESCA_DE.txt", 
                          "GBM" = "data/TCGA_GBM_DE.txt", 
                          "HNSC" = "data/TCGA_HNSC_DE.txt", 
                          "KICH" = "data/TCGA_KICH_DE.txt", 
                          "KIRC" = "data/TCGA_KIRC_DE.txt", 
                          "KIRP" = "data/TCGA_KIRP_DE.txt", 
                          "LIHC" = "data/TCGA_LIHC_DE.txt", 
                          "LUAD" = "data/TCGA_LUAD_DE.txt", 
                          "LUSC" = "data/TCGA_LUSC_DE.txt", 
                          "PCPG" = "data/TCGA_PCPG_DE.txt", 
                          "PRAD" = "data/TCGA_PRAD_DE.txt", 
                          "READ" = "data/TCGA_READ_DE.txt", 
                          "STAD" = "data/TCGA_STAD_DE.txt", 
                          "THCA" = "data/TCGA_THCA_DE.txt", 
                          "UCEC" = "data/TCGA_UCEC_DE.txt",
                          "GBM.sig" = "data/TCGA_GBM_Signature.txt",
                          "CRC.sig" = "data/TCGA_CRC_Signature.txt",
                          "BRCA.sig" = "data/TCGA_BRCA_Signature.txt",
                          "PDX.GBM1" = "data/Table_G1_1_PDX_Group1_L1000_only.txt",
                          "PDX.GBM2" = "data/Table_G2_1_PDX_Group2_L1000_only.txt",
                          "PDX.GBM3" = "data/Table_G3_1_PDX_Group3_L1000_only.txt" ,
                          "PDX.GBM4" = "data/Table_G4_1_PDX_Group4_L1000_only.txt") 
  
  if (is.null(options_Bioactiv)==TRUE){
    input_bioactivities <- "IC50"
    message('No drug activity measurement specified. Defaulting to IC50. ')
  } else {
    
    input_bioactivities <- options_Bioactiv
  }
  
  input_bioactivFile <- switch(input_bioactivities,
                                "IC50" = "data/Bioactivities_IC50.txt",
                                "Kd" = "data/Bioactivities_Kd.txt",
                                "Ki" = "data/Bioactivities_Ki.txt",
                                "Potency" = "data/Bioactivities_Potency.txt"
  )
  
  Bioactivities <- read.table(file=input_bioactivFile, sep ="\t", header=TRUE)
  
  TCGA_Sig <- read.csv(file=input_disease, sep="")
  
  JQ1_Sig_v <- JQ1_Sig
  JQ1_Sig_v <- data.frame(JQ1_Sig_v, Genes=row.names(JQ1_Sig_v))
  colnames(JQ1_Sig_v) <- c("JQ1","Genes")
  TCGA_Sig_v <- TCGA_Sig
  V <- unique(as.numeric(as.character(JQ1_Sig_v$JQ1))) 
  V2 <- max(abs(V))
  V4 <- as.numeric(n_bins)/100
  V5 <- V2*V4 
  V6 <- round(V5)
  JQ1_Sig_v$JQ1 <- as.numeric(as.character(JQ1_Sig_v$JQ1))
  JQ1_Sig_v2 <- JQ1_Sig_v[which(JQ1_Sig_v$JQ1 > V6 | JQ1_Sig_v$JQ1 < -V6),]
  JQ1_Sig_v2 <- JQ1_Sig_v2[which(JQ1_Sig_v2$JQ1 > 0 | JQ1_Sig_v2$JQ1 < 0),]
  jq1_genes <- as.character(JQ1_Sig_v2$Genes)
  Drugs_Sigs2 <- Drugs_Sigs[,jq1_genes]
  tt <-  apply(Drugs_Sigs2, 1, function(x){x*JQ1_Sig_v2$JQ1})
  
  
  # calculate similarity of drug signatures to selected reference drug (ratio of genes that are concordant with the reference drug signature to number of genes that are discordant to reference drug signature)
  tt2 <-  apply(tt, 2, function(x) { 
  a <- sum(x>0)
  b <- sum(x<0)
  b[b==0] <- 1 # replace b with 1 so we can divide by that number. Another way would be to add one to both a and b
  c <- a/b
  #c
  return(c(a,b,c))
  })
  tt3 <- as.data.frame(t(tt2))
 
  ### ### ### ### ### ### ### ### ### ### ### ###
  # calculate ratio of genes that are discordant to the Disease Signature (and are not affected by the selected reference drug) to genes that are concordant with the Disease Signature (and non-reference drug signatures)
  
  jq1_genes2 <- as.character(JQ1_Sig_v[which(JQ1_Sig_v$JQ1==0),2])
  
  genes <- colnames(Drugs_Sigs)
  TCGA_Sig_v$Genes.log2FoldChangelog2FoldChange <- as.numeric(TCGA_Sig_v$log2FoldChange)
  TCGA_Sig_v2 <- TCGA_Sig_v[,-1,drop=FALSE]
  TCGA_Sig_v3 <- aggregate(TCGA_Sig_v2, by = list(TCGA_Sig_v$Genes),FUN=mean) # this is in case we have duplicate gene symbols in the signature
  non_jq1_genes <- intersect(as.character(TCGA_Sig_v3$Group.1),jq1_genes2)
  Drugs_Sigs3 <- Drugs_Sigs[,non_jq1_genes]
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
  Final <- merge(pp3[,3,drop=F],tt3[,3,drop=F],by="row.names")
  
  colnames(Final) <- c("Drug","Disease_Discordance", "Reference_Drug_Orthogonality")
  
  Final
  Final[,2] <- round(Final[,2],digits=3)
  Final[,3] <- round(Final[,3],digits=3)
  
  bioact.moa <- merge(Bioactivities, SM_MOA, by = "Drugs")
  f3 <- merge(Final, bioact.moa, by.x = "Drug", by.y = "Drugs", all.x = TRUE)
  
  f3.ord <-f3[with(f3, order(-Disease_Discordance, Reference_Drug_Orthogonality)), ]
  
  if (is.null(molecTarget)==TRUE){
    target <- as.character(f3.ord$target_gene_symbol[1])
    message(paste0('No molecular target of interest specified. Defaulting to ', target, '. '))
  } else {
    target <- molecTarget
  }
  
  drugs.sameTargets1 <- na.omit(as.data.frame(f3.ord[f3.ord$target_gene_symbol==target,]))
  
  colnames(drugs.sameTargets1)[1] <- "Drug"
  
  # library(RColorBrewer)
  
  pal <- brewer.pal(nrow(drugs.sameTargets1),"BrBG")
  
  SM_MOA <- read.csv(file="data/L1000_SM_MOA.csv", header=TRUE)
  
  Final2 <- merge(Final, SM_MOA, by.x = "Drug", by.y = "Drugs", all.x = TRUE)
  
  p1 <- plot_ly(Final, x = ~Reference_Drug_Orthogonality, y = ~Disease_Discordance, type = 'scatter', alpha=0.7, marker = list(size = 14),
                mode = 'markers', hoverinfo = 'text', text = ~paste(Drug,'<br>', "Ratio:", Disease_Discordance, Reference_Drug_Orthogonality)) %>% layout(dragmode = "select")
  
  plotly.plot <- p1 %>%
    add_data(drugs.sameTargets1) %>%
    
    add_markers(x = ~Reference_Drug_Orthogonality, y = ~Disease_Discordance, color = ~Median, colors = pal) %>% colorbar(title = paste0("Median ", input_bioactivities))
  
  res2 <- list(assign("res.table", Final2), assign("res.plotly", plotly.plot))
  
  return(res2)
}


#' Drugs_SigsR.class
#'
#' Creates a Drugs_SigsR-class object from a list containing a drug synergy table and a plotly plot.
#' Returns a Drugs_SigsR-class object. 
#'
#' @keywords internal
#'
#' @export


