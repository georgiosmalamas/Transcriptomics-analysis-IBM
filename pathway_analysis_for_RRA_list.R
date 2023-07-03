##libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("RobustRankAggreg")
install.packages("stringr")
install.packages("visNetwork")
install.packages("xml2")
install.packages("igraph")
install.packages("forcats")
install.packages("scales")
install.packages("ggplot2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("biomaRt")
BiocManager::install("GO.db")
BiocManager::install("hmdbQuery")
BiocManager::install("enrichplot")
library(igraph)
library(visNetwork) 
library("RobustRankAggreg")
library("httr")
library("clusterProfiler")
library("org.Hs.eg.db")
library("biomaRt")
library("GO.db")
library("stringr")
library("enrichplot")
library("xml2")
library("hmdbQuery")
library("ggplot2")
library("scales")
library("forcats")
# RRA

## load all lists from excel
df1_DE_Genes <- read_excel("df1_DE_Genes.xlsx")
df2_DE_Genes <- read_excel("df2_DE_Genes.xlsx")
df3_DE_Genes <- read_excel("df3_DE_Genes.xlsx")
df4_DE_Genes <- read_excel("df4_DE_Genes.xlsx")
df5_DE_Genes <- read_excel("df57_DE_Genes.xlsx")

## Order the genes based on adj.p-value to prepare for RRA
ordf1 <- df1_DE_Genes[order(df1_DE_Genes$adj.P.Val),]
ordf2 <- df2_DE_Genes[order(df2_DE_Genes$adj.P.Val),]
ordf3 <- df3_DE_Genes[order(df3_DE_Genes$adj.P.Val),]
ordf4 <- df4_DE_Genes[order(df4_DE_Genes$adj.P.Val),]
ordf5 <- df5_DE_Genes[order(df5_DE_Genes$padj),]

## Make sample input data only for the gene symbols
glist <- list(as.vector(ordf4$hgnc_symbol), as.vector(ordf1$SYMBOL), as.vector(ordf3$hgnc_symbol) , as.vector(ordf2$hgnc_symbol), as.vector(ordf5$genes))

## Get rid of the NAs
for (i in 1:length(glist)) {
  glist[[i]] = na.omit(glist[[i]])
}

## Robust Rank Aggregation analysis

aR <- aggregateRanks(glist=glist)
colnames(aR) <- c("Gene_Symbol", "RRA_scores")

### FDR correction

aR$adj.pval <- aR$RRA_scores*length(glist)
ranks_result <- aR[aR$adj.pval <=0.05,]

# GO analysis

rG <-as.data.frame(ranks_result$Gene_Symbol)
gene <- rG[,1]
gene_conversion <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego <- enrichGO(gene_conversion$ENTREZID, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, readable=T)
res1GO <-ego@result
ego2 <- enrichGO(gene_conversion$ENTREZID, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="MF", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, readable=T)
res2GO <-ego2@result
ego3 <- enrichGO(gene_conversion$ENTREZID, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="CC", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, readable=T)
res3GO <-ego3@result

resadjpf1 <- subset(res1GO, p.adjust <0.05)
resadjpf2 <- subset(res2GO, p.adjust <0.05)
resadjpf3 <- subset(res3GO, p.adjust <0.05)

### Create txt with the enrichment resutls
write.table(resadjpf$ID, "GOtermsBP.txt",row.names=FALSE, quote = FALSE)
write.table(resadjpf$ID, "GOtermsMF.txt",row.names=FALSE, quote = FALSE)
write.table(resadjpf$ID, "GOtermsCC.txt",row.names=FALSE, quote = FALSE)

## Wikipathways

WP <- enrichWP(gene_conversion$ENTREZID, organism = "Homo sapiens")
resWP <- DOSE::setReadable(WP, org.Hs.eg.db, keyType = "ENTREZID")
resWP1 <-resWP@result
write.table(resWP1$ID, "WPterms.txt",row.names=FALSE, quote = FALSE)
barplot(resWP, showCategory = 20)
dotplot(resWP, showCategory = 20)

## Find mRNA terms for WP and GO(It can be done for any terms) 
### For WP
query_wikipathways <- function(query) {
  url <- "https://webservice.wikipathways.org/findPathwaysByLiterature"
  params <- list(query = query, format = "json")
  res <- httr::GET(url, query = params)
  httr::content(res, as = "parsed", type = "application/json")$result
}

query <- "mRNA splicing"
results <- query_wikipathways(query)
pw_names <- sapply(results, function(x) x$name)
pw_urls <- sapply(results, function(x) x$url)
WPtermsmRNA <- data.frame(Name = pw_names, URL = pw_urls)
WPtermsmRNA$WPterms <-str_split_fixed(WPtermsmRNA$URL, ":", 3)[,3]
intersect(WPtermsmRNA$WPterms,resWP1$ID)
mRNAsplWP <- resWP1[intersect(WPtermsmRNA$WPterms,resWP1$ID),]
GO <- as.list(GOTERM)
select <-which(grepl('mRNA splicing', GO))

GO_id = NULL
for (i in 1:length(select)) {
  GO_id <- cbind (GO_id, GO[[select[i]]]@GOID)
  
}

GO_term = NULL
for (i in 1:length(select)) {
  GO_term <- cbind (GO_term, GO[[select[i]]]@Term)
  
}
GO_ont = NULL
for (i in 1:length(select)) {
  GO_ont <- cbind (GO_ont, GO[[select[i]]]@Ontology)
  
}
GO_def = NULL
for (i in 1:length(select)) {
  GO_def <- cbind (GO_def, GO[[select[i]]]@Definition)
  
}

mRNAsplicing_GOterms <- data.frame(t(GO_id),t(GO_term),t(GO_def),t(GO_ont))
colnames(mRNAsplicing_GOterms) <- c("ID","Term","Definition","Ontology")
intersect(mRNAsplicing_GOterms$ID,rownames(resallGO))

##  Visualization of the enrtichment results
dotplot(ego, showCategory=30) + ggtitle("dotplot for GO")
dotplot(WP, showCategory=30) + ggtitle("dotplot for WP")


# PPI

## Construct the URL for the STRING API request
string_api_get_enrichment <- function(api_base_url, protein_list, species=9606, caller_identity=NULL) {
  
  url <- sprintf("%s/enrichment?identifiers=%s&species=%d", api_base_url, paste(protein_list, collapse="%0d"), species)
  if (!is.null(caller_identity)) {
    url <- sprintf("%s&caller_identity=%s", url, caller_identity)
  }
  
  ## Send the request to the STRING API and parse the response
  response <- httr::GET(url)
  if (response$status_code != 200) {
    stop("Error: Could not retrieve enrichment data from STRING API.")
  }
  results <- jsonlite::fromJSON(httr::content(response, as="text"))
  
  ## Extract the Ensembl IDs from the API response
  ensembl_ids <- sapply(results$terms, function(x) {
    unlist(strsplit(x$id, "\\."))[1]
  })
  
  return(ensembl_ids)
}


### Plot the interactions

#### nodes and edges
nodes <- unique(c(interactions$Protein1, interactions$Protein2))
nodes_df <- data.frame(id = nodes, label = nodes)

edges_df <- data.frame(from = interactions$Protein1, to = interactions$Protein2, width = interactions$Prob.numbers)

#### Create the interactive graph
visNetwork(nodes_df, edges_df) %>%
  visEdges(arrows = "to") %>%
  visPhysics(solver = "forceAtlas2Based") %>%
  visOptions(highlightNearest = TRUE, selectedBy = "label")

### Otherwise you can use the igraph package
g <- igraph::graph_from_data_frame(interactions[,c('Protein1','Protein2')])
E(g)$weight <- interactions$Prob.numbers
V(g)$type <- ifelse(V(g)$name %in% E(g))
igraph::tkplot(g)

## Metabolomics integration 
###First find KEGG ids


library(MetaboSignal)
keggID <- MS_convertGene(gene_conversion$SYMBOL, organism_name = "human", organism_code = "hsa")

#### Retrieve the HMDB metabolites from your gene list

GeneNames_HMDB_metabolites <- read.csv("GeneNames.csv")
GeneNames_HMDB_metabolites <- na.omit(GeneNames_HMDB_metabolites)
hmdb_id <- GeneNames_HMDB_metabolites$HMDP_ID
hmdb_metabolites_for_each_gene = list()
### Take for each gene name
for (i in 1:length(hmdb_id)) {
  hmdb_metabolites_for_each_gene[[i]] <- HmdbEntry(prefix = "http://www.hmdb.ca/proteins/",
                                                   id = hmdb_id[i],
                                                   keepFull = TRUE)
  
}

hl1 <- hmdb_metabolites_for_each_gene[[1]]@store$metabolite_associations
hl2 <- hmdb_metabolites_for_each_gene[[2]]@store$metabolite_associations
hl3 <- hmdb_metabolites_for_each_gene[[3]]@store$metabolite_associations
hl4 <- hmdb_metabolites_for_each_gene[[4]]@store$metabolite_associations
hl5 <- hmdb_metabolites_for_each_gene[[5]]@store$metabolite_associations
hl6 <- hmdb_metabolites_for_each_gene[[6]]@store$metabolite_associations
hl7 <- hmdb_metabolites_for_each_gene[[7]]@store$metabolite_associations
hl8 <- hmdb_metabolites_for_each_gene[[8]]@store$metabolite_associations
hl9 <- hmdb_metabolites_for_each_gene[[9]]@store$metabolite_associations
hl10 <- hmdb_metabolites_for_each_gene[[10]]@store$metabolite_associations
hl11 <- hmdb_metabolites_for_each_gene[[11]]@store$metabolite_associations
hl12 <- hmdb_metabolites_for_each_gene[[12]]@store$metabolite_associations
hl13 <- hmdb_metabolites_for_each_gene[[13]]@store$metabolite_associations
hl14 <- hmdb_metabolites_for_each_gene[[14]]@store$metabolite_associations
hl15 <- hmdb_metabolites_for_each_gene[[15]]@store$metabolite_associations
hl16 <- hmdb_metabolites_for_each_gene[[16]]@store$metabolite_associations
hl17 <- hmdb_metabolites_for_each_gene[[17]]@store$metabolite_associations
hl18 <- hmdb_metabolites_for_each_gene[[18]]@store$metabolite_associations
### Make each of it as dataframe
hl1test <- as.data.frame(unlist(hl1[2,]))
hl2test <- as.data.frame(unlist(hl2[2,]))
hl3test <- as.data.frame(unlist(hl3[2,]))
hl4test <- as.data.frame(unlist(hl4[2,]))
hl5test <- as.data.frame(unlist(hl5[2,]))
hl6test <- as.data.frame(unlist(hl6[2,]))
hl7test <- as.data.frame(unlist(hl7[2,]))
hl8test <- as.data.frame(unlist(hl8[2,]))
hl9test <- as.data.frame(unlist(hl9[2,]))
hl10test <- as.data.frame(unlist(hl10[2,]))
hl11test <- as.data.frame(unlist(hl11[2,]))
hl12test <- as.data.frame(unlist(hl12[2,]))
hl13test <- as.data.frame(unlist(hl13[2,]))
hl14test <- as.data.frame(unlist(hl14[2,]))
hl15test <- as.data.frame(unlist(hl15[2,]))
hl16test <- as.data.frame(unlist(hl16[2,]))
hl17test <- as.data.frame(unlist(hl17[2,]))
hl18test <- as.data.frame(unlist(hl18[2,]))


colnames(hl1test) <- "Metabolites"
colnames(hl2test) <- "Metabolites"
colnames(hl3test) <- "Metabolites"
colnames(hl4test) <- "Metabolites"
colnames(hl5test) <- "Metabolites"
colnames(hl6test) <- "Metabolites"
colnames(hl7test) <- "Metabolites"
colnames(hl8test) <- "Metabolites"
colnames(hl9test) <- "Metabolites"
colnames(hl10test) <- "Metabolites"
colnames(hl11test) <- "Metabolites"
colnames(hl12test) <- "Metabolites"
colnames(hl13test) <- "Metabolites"
colnames(hl14test) <- "Metabolites"
colnames(hl15test) <- "Metabolites"
colnames(hl16test) <- "Metabolites"
colnames(hl17test) <- "Metabolites"
colnames(hl18test) <- "Metabolites"
### Bind them together
metabolm <- rbind(hl18test,hl17test,hl16test,hl15test,hl14test,hl13test,hl12test,hl11test,hl10test,hl9test,hl8test,hl7test,hl6test,hl5test,hl4test,hl3test,hl2test,hl1test)
metabolmuniq <- unique(metabolm) ##unique metabolites for the MetaboAnlyst tool



###circos plot

##ChordDiagramm
newgen_metab <- data.frame (GeneNames=GeneNames_HMDB_metabolites$GeneName,Metabolites)
metabolitesnew = strsplit(newgen_metab$Metabolites,",")
df2 = data.frame(
  gene = rep(newgen_metab$GeneNames, sapply(metabolitesnew, length)),
  metabolite = unlist(metabolitesnew),
  stringsAsFactors = FALSE
)
df3 <-df2[,1:2]
col_order <- c("metabolite", "gene")
df4 <-df3[,col_order]
circos.info()
circos.clear()
df5 <- df4[c(1:29,68),]
df6 <- df3[c(1:29,68),]


chordDiagram(df6, annotationTrack = c("gene"), 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df6))))),list(bg.border = "purple"))



## Skeletal muscle expression from Human Atlas platform

### Select for muscle and for the  genes(first load from rnagtex file on the playform itself)

newdata <- subset(rna_tissue_gtex, Tissue== 'skeletal muscle')
newdf <- newdata[newdata$`Gene name` %in% GeneNames$X1,]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

### Set the order and subset the genes
k <-order(newdf$pTPM)
newdf <- newdf[k,]
dfreor <- forcats::fct_reorder(newdf$`Gene name`, -newdf$pTPM)
newdfmet <-newdf[newdf$'Gene name' %in% genes_metabolites$GeneNames,]

### Create a histogram with bars 
ggplot(newdf, aes(x = forcats::fct_reorder(newdf$`Gene name`, newdf$pTPM), y = newdf$pTPM, fill= newdf$'Gene name')) +
  geom_col(fill= rainbow(74), position= "dodge") +
  facet_wrap(~Tissue) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(title = "Transcript expression level on muscle tissue",
       x = "Genes",
       y = "pTPM")

### only for metabolites

ggplot(newdfmet , aes(x = forcats::fct_reorder(newdfmet$`Gene name`, newdfmet$pTPM), y = newdfmet$pTPM, fill= newdfmet$'Gene name')) +
  geom_col(fill= rainbow(18), position= "dodge") +
  facet_wrap(~Tissue) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(title = "Transcript expression level on muscle tissue for metabolites-associated genes",
       x = "Genes",
       y = "pTPM")

### Explore the ratio

newdfratio <- rna_tissue_gtex[rna_tissue_gtex$`Gene name` %in% GeneNames$X1,]
newdf$alltissuepTPM <- NA
result <- aggregate(pTPM ~ 'Gene name', newdfratio, sum)

for (i in 1:74) {
    lop<- grep(newdf$'Gene name'[i], newdfratio$`Gene name`)
    newdf$alltissuepTPM[i] <-sum(newdfratio$pTPM[lop])
}

newdf$muscleratiopTPM <- newdf$pTPM*100/newdf$alltissuepTPM  
ggplot(newdf, aes(x = forcats::fct_reorder(newdf$`Gene name`, newdf$muscleratiopTPM), y = newdf$muscleratiopTPM, fill= newdf$'Gene name')) +
  geom_col(fill= rainbow(74), position= "dodge") +
  facet_wrap(~Tissue) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(title = "Transcript expression ratio",
       x = "Genes",
       y = "pTPM %")  

### only for metabolites-associated genes 
ggplot(newdfmet , aes(x = forcats::fct_reorder(newdfmet$`Gene name`, newdfmet$muscleratiopTPM), y = newdfmet$muscleratiopTPM, fill= newdfmet$'Gene name')) +
  geom_col(fill= rainbow(18), position= "dodge") +
  facet_wrap(~Tissue) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(title = "Transcript expression ratio for metabolites-associated genes",
       x = "Genes",
       y = "pTPM ratio")

