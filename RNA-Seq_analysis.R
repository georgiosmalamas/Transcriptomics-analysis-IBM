##Libraries

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("BiocStyle")
BiocManager::install("limma")
BiocManager::install("oligo")
BiocManager::install("biomaRt")
BiocManager::install("DESeq2")
install.packages("openxlsx")
install.packages("data.table")
install.packages("stringr")

library(BiocStyle)
library(dplyr)
library(tidyr)
library(limma)
library(oligo)
library(openxlsx)
library(GEOquery)
library(DESeq2)
library(data.table)
library(stringr)
library(biomaRt)



data <- "GSE151757"
cohort_control_name <- "Amputee"
cohort_case_name <- "Inclusion body myositis"

perform_de_analysis <-function (data,cohort_case_name,cohort_control_name){
  
  if (file.exists(data)== TRUE) {
  df<-read.delim(data) #Formal expression set
  datainfo <- Biobase::pData(df) #datainfo
  columnname <- datainfo[, c(which(grepl(cohort_control_name , datainfo)))]
  datainfo <- subset(datainfo, columnname[ncol(columnname)]== cohort_case_name|columnname[ncol(columnname)]== cohort_control_name) #datainfo for disease 
  ##finding the normalized counts
  files <- list.files(path = ".", full.names = TRUE)
  selected_files <- files[ grepl("normalized", basename(files))]
  b2 <-read.csv(selected_files,sep = "")
  k = datainfo[,1]
  dat <- b2[,k] 
  extra_col <- as.data.frame(str_split_fixed(row.names(dat), ";",2))
  
  } else if (file.exists(data)!= TRUE) {
  df <- getGEO(data)[[1]]
  datainfo <- Biobase::pData(df)
  columnname <- datainfo[, c(which(grepl(cohort_control_name , datainfo)))]
  datainfo <- subset(datainfo, columnname[ncol(columnname)]== cohort_case_name|columnname[ncol(columnname)]== cohort_control_name)  
  
  sfiles = getGEOSuppFiles(data)
  fnames = rownames(sfiles)
  # take the supplemental file for normalized counts
  b2 = read.csv(fnames[1],sep="")
  k = datainfo[,1]
  dat <- b2[,k] 
  
  extra_col <- as.data.frame(str_split_fixed(row.names(dat), ";",2))
  
  
  } else {
    print("No analysis for you")
  }
  
  ## Creating metaData object 
  md1 <- which(grepl(cohort_case_name , datainfo))[2]
  md2 <- which(grepl("GSM" , datainfo))
  
  metadata <- datainfo[,c(md1,md2)]
  colnames(dat) <-metadata[,ncol(metadata)]
  
  design <- as.character(metadata[,1])
  colnames(metadata)[1] <- "design"
  ## The analysis
  
  dds <- DESeqDataSetFromMatrix(countData = round(dat), colData = metadata, design = ~design)
  dds <- DESeq(dds)
  res <- results(dds)
  res2 <-as.data.frame(res)
  
  
  ## Selecting protein coding RNAs
  
  ensembl <- useMart("ensembl")
  ensembl = useDataset(dataset="hsapiens_gene_ensembl",mart=ensembl)
  
  lookuptable <- getBM(
    mart = ensembl,
    attributes = c(
      'ensembl_gene_id',
      'hgnc_symbol',
      'gene_biotype'),
    uniqueRows = TRUE)
  
  protein_coding_ENSG <- lookuptable$ensembl_gene_id[lookuptable$gene_biotype == "protein_coding"]
  
  valid <- res2[!is.na(res2$padj), ]
  colvalid <- substr(rownames(valid), 1, 15)
  valid <- valid[substr(rownames(valid), 1, 15) %in% protein_coding_ENSG, ]
  
  dysregulated <- valid[valid$padj < 0.01, ]
  upregulated <- dysregulated[dysregulated$log2FoldChange >= 1.5, ]
  downregulated <- dysregulated[dysregulated$log2FoldChange <= -1.5, ]
  DEG <- rbind(upregulated,downregulated)
  w= write.xlsx(DEG, 'DEGs.xlsx')
  return(w)
}


perform_de_analysis("GSE151757","Inclusion body myositis","Amputee")