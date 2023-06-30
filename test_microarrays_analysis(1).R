## LIbraries

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("BiocStyle")
BiocManager::install("limma")
BiocManager::install("oligo")
install.packages("openxlsx")

library(GEOquery)
library(BiocStyle)
library(dplyr)
library(tidyr)
library(limma)
library(oligo)
library(openxlsx)


## Function

perform_de_analysis <-function (data,cohort_case_name,cohort_control_name) {
  
  
  
  if (file.exists(data)== TRUE) {
    SDRF <- read.delim(data)
    
    column_data<-SDRF[, c(which(grepl(cohort_case_name, SDRF)))]
    
    SDRF <-  subset(SDRF, column_data[ncol(column_data)]== cohort_case_name|column_data[ncol(column_data)]== cohort_control_name)
    
    columnnames <- SDRF[, c(which(grepl("CEL" , SDRF)))]
    
    SDRF <- AnnotatedDataFrame(SDRF)
    
    
    raw_data <- oligo::read.celfiles(columnnames[1], verbose = FALSE, phenoData = SDRF)
    
    stopifnot(validObject(raw_data))
    
    palmieri_eset_norm <- oligo::rma(raw_data, target = "core")
    
    
    design <- as.character(Biobase::pData(palmieri_eset_norm)[,ncol(SDRF)])
    design_palmieri <-model.matrix(~0 + design)
    colnames(design_palmieri) = c("disease","control")
    
    
    fit <- lmFit(palmieri_eset_norm,design_palmieri)
    
  } else if (file.exists(data)!= TRUE) {
    
    expdf <-getGEO(data)[[1]]
    SDRF <- Biobase::pData(expdf)
    
    columnname <- SDRF[, c(which(grepl(cohort_case_name , SDRF)))]
    
    SDRF <- subset(SDRF, columnname[ncol(columnname)]== cohort_case_name|columnname[ncol(columnname)]== cohort_control_name)  
    
    
    expr <-Biobase::exprs(expdf)
    expr_disease <- expr[,c(row.names(SDRF))]
    
    columnname2 <-which(grepl(cohort_case_name , SDRF))[2]
    colnames(expr_disease) = t(SDRF[columnname2])
    
    design<- as.character(colnames(expr_disease))
    design_palmieri <-model.matrix(~0 + design)
    colnames(design_palmieri) = c("control","disease")
    
    
    fit <- lmFit(expr_disease,design_palmieri)
    
  } else {
    print("No analysis for you")
  }
  
  
  
  
  
  contrast_matrix <- makeContrasts( disease - control, levels = design_palmieri)
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  table_diseaseDEG <- topTable(fit2, adjust.method="BH", n=Inf, p.value=0.05)
  
  table_diseaseDEG$namesID =rownames(table_diseaseDEG)
  
  w= write.xlsx(table_diseaseDEG, 'DEGs.xlsx')
  return(w)
  
}

##Run the function

perform_de_analysis(data= "GSE48280",cohort_case_name= "inclusion body myositis",cohort_control_name="healthy control") # GEO example