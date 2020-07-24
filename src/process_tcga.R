## This is script for processing TCGA data. 
# a) select studies with greater than 30 samples (all cancer types have more than 30 samples)
# b) split expression, gene and other metadata into different Rds objects

library(recount)
library(dplyr)

dNames <- commandArgs(TRUE)
print(dNames)
inputDir <- dNames[1]
outputDir <- dNames[2]

# Read in TCGA data
tcga <- readRDS(paste0(inputDir, "/tcga.Rds"))
expr.tcga <- as.data.frame(tcga@assays@data$counts)
sample_metadata.tcga <- data.frame(colData(tcga))
gene_data.tcga <- data.frame(rowData(tcga))
rownames(gene_data.tcga) <- c(1:dim(gene_data.tcga)[1])

sample.size.tissue <- table(tcga$gdc_cases.project.project_id)
all(sample.size.tissue >= 30)
gene_id <- c()
bp_length <- c()
gene_symbol <- c()
for(i in c(1:dim(gene_data.tcga)[1])){
  gene_id[i] <- gene_data.tcga[i, 1]
  bp_length[i] <- gene_data.tcga[i, 2]
  gene_symbol[i] <- unlist(gene_data.tcga[i, 3])[1]
}

gene_data.tcga <- data.frame(gene_id, bp_length, gene_symbol)
gene_data.tcga <- distinct(gene_data.tcga, gene_symbol, .keep_all = TRUE)
gene_data.tcga <- gene_data.tcga[!is.na(gene_data.tcga$gene_symbol), ]

# Subset the expression data
expr.tcga <- expr.tcga[rownames(expr.tcga) %in% gene_data.tcga$gene_id, ]
expr.tcga <- expr.tcga[match(gene_data.tcga$gene_id, rownames(expr.tcga)), ]
rownames(expr.tcga) <- gene_data.tcga$gene_symbol

saveRDS(expr.tcga, paste0(outputDir, "/tcga/expr_tcga.rds"))
saveRDS(sample_metadata.tcga, paste0(outputDir, "/tcga/sample_metadata.rds"))
saveRDS(gene_data.tcga, paste0(outputDir, "/tcga/gene_data.rds"))
