rm(list= ls())

#dirName <- "/work-zfs/abattle4/prashanthi/consensus_networks/data/sra_metadata/"
#rawExp <- "/work-zfs/abattle4/prashanthi/consensus_networks/data/rpkm/sra.Rds"
#saveDir <- "/work-zfs/abattle4/prashanthi/consensus_networks/data/rpkm/"
#datDir <- "/work-zfs/abattle4/prashanthi/consensus_networks/data/"
#dir.create(dirName)
################################################################################################################################
## Get inputargs
inputArgs <- commandArgs(TRUE)
dirName <- inputArgs[1]
rawExp <- inputArgs[2]
datDir <- inputArgs[3]
saveDir <- inputArgs[4]

## Create directories
dir.create(dirName)
## processed recountDir
processed_recountDir <- paste0(datDir, "/expr_data/recount/") 
dir.create(processed_recountDir)
## processed gtexDir
processed_gtexDir <- paste0(datDir, "/expr_data/GTEx/")
dir.create(processed_gtexDir)

################################################################################################################################

library(SRAdb)
library(dplyr)
library(plyr)
library(recount)
sra <- readRDS(rawExp)
runids <- sra$run # get run ids

# get corresponding metadata from SRAdb
dbFilename <- paste(dirName, "SRAmetadb.sqlite", sep = "")
#sqlfile <- getSRAdbFile(, destfile = paste(dbFilename, ".gz", sep = "")) ## not working on marcc right now, using the previously downloaded version of db
sra_con <- dbConnect(SQLite(),dbFilename)
res <- dbGetQuery(sra_con, paste("select * from sra where run_accession in ('",paste(runids, collapse = "','"),"')", sep = ""))

run_diff <- setdiff(runids, res$run_accession)
study_diff <- table(sra$project[sra$run %in% run_diff])

## only select runs that have annotations/metadata
sra <- sra[,match(res$run_accession,sra$run)]

potential_cell_lines <- c(grep("cell_line", res$sample_attribute), grep("cells", res$experiment_title))
# res$experiment_title[potential_cell_lines]

smallrna_samples <- res[res$library_selection %in% "size fractionation",]$run_accession
sra <- sra[,!sra$run %in% smallrna_samples]

studies <- unique(sra$project)
# replicate_exp <- table(sra$experiment)>1
# studies <- unique(sra$project[sra$experiment %in% names(replicate_exp)[which(replicate_exp)]])

each_study_count <- lapply(studies, function(x,y){
  print(x)
  this_y <- y[, y$project %in% x]
  this_study <- as.data.frame(t(this_y@assays@data$counts))
  this_study$exp <- this_y$experiment
  if(length(unique(this_y$experiment)) == length(unique(this_y$run))){
    this_count <- this_study
  }else{
    this_count <- this_study %>% group_by(exp) %>% summarise_all(median) %>% as.data.frame
  }
  this_count
}, sra)

each_study_count <- lapply(each_study_count, function(this_study){
  rownames(this_study) <- this_study$exp
  this_study$exp <- NULL
  this_study
})

names(each_study_count) <- studies


## Exclude samples with more than 50% of genes with 0 expression value
total_genes = ncol(each_study_count[[1]])

prop_genes_notexp = lapply(each_study_count, function(x) apply(x, 1, function(y) (sum(y<=0)/total_genes)<0.5))

each_study_count = mapply(function(x,y){
  x[y,]
}, each_study_count, prop_genes_notexp, SIMPLIFY = FALSE)

saveRDS(each_study_count, file = paste(datDir,"/replicates_merged_sra/rpkm_replicates_merged.Rds", sep = ""))

## Merge the replicates to get one single matrix 
expr.recount <- ldply(each_study_count, data.frame)
expr.recount$.id <- NULL
exp.id <- c()
for(iexp in each_study_count){
  exp.id <- c(exp.id, rownames(iexp))
}
row.names(expr.recount) <- exp.id
expr.recount <- t(expr.recount)
saveRDS(expr.recount, paste0(saveDir,"expr_recount.rds"))

## Save gene and sample metadata 
sra <- add_predictions(sra)
sample_metadata <- data.frame(colData(sra))
sample_metadata <- sample_metadata[sample_metadata$experiment %in% colnames(expr.recount), ]
sample_metadata <- distinct(sample_metadata, experiment, .keep_all = TRUE)
sample_metadata <- sample_metadata[match(colnames(expr.recount), sample_metadata$experiment), ]
gene_data <- data.frame(rowData(sra))
rownames(gene_data) <- c(1:dim(gene_data)[1])
saveRDS(sample_metadata, paste0(saveDir, "sample_metadata.rds"))
saveRDS(gene_data, paste0(saveDir, "gene_data.rds"))

# Read in projects to be excluded 
projects.excl <- readRDS(paste0(datDir, "projects_excluded.rds"))
sample_metadata <- sample_metadata[!sample_metadata$project %in% projects.excl$accession, ]
sample.size <- table(sample_metadata$project)
small_projects <- names(sample.size)[sample.size < 30]

sample_metadata <- sample_metadata[!sample_metadata$project %in% small_projects, ]

# Subset only to these samples and studies 
expr.recount <- expr.recount[ ,colnames(expr.recount) %in% sample_metadata$experiment]

# For the gene data we will select the first symbol 
gene_id <- c()
bp_length <- c()
gene_symbol <- c()
for(i in c(1:dim(gene_data)[1])){
  gene_id[i] <- gene_data[i, 1]
  bp_length[i] <- gene_data[i, 2]
  gene_symbol[i] <- unlist(gene_data[i, 3])[1]
}

gene_data <- data.frame(gene_id, bp_length, gene_symbol)
gene_data <- distinct(gene_data, gene_symbol, .keep_all = TRUE)

# Subset the expression data
expr.recount <- expr.recount[rownames(expr.recount) %in% gene_data$gene_id, ]
expr.recount <- expr.recount[match(gene_data$gene_id, rownames(expr.recount)), ]
rownames(expr.recount) <- gene_data$gene_symbol


saveRDS(expr.recount, paste0(processed_recountDir,"/expr_recount.rds"))
saveRDS(sample_metadata, paste0(processed_recountDir,"/sample_metadata.rds"))
saveRDS(gene_data, paste0(processed_recountDir,"/gene_data.rds"))

# Read in GTEx data
GTEx <- readRDS(paste0(datDir,"/rpkm/SRP012682.Rds"))
expr.gtex <- as.data.frame(GTEx@assays@data$counts)
sample_metadata.gtex <- data.frame(colData(GTEx))
gene_data.gtex <- data.frame(rowData(GTEx))
rownames(gene_data.gtex) <- c(1:dim(gene_data.gtex)[1])

sample.size.tissue <- table(GTEx$smtsd)
select.tissues <- names(sample.size.tissue)[sample.size.tissue >= 30]
sample_metadata.gtex <- sample_metadata.gtex[sample_metadata.gtex$smtsd %in% select.tissues, ]
expr.gtex <- expr.gtex[ ,sample_metadata.gtex$run]

gene_id <- c()
bp_length <- c()
gene_symbol <- c()
for(i in c(1:dim(gene_data.gtex)[1])){
  gene_id[i] <- gene_data.gtex[i, 1]
  bp_length[i] <- gene_data.gtex[i, 2]
  gene_symbol[i] <- unlist(gene_data.gtex[i, 3])[1]
}

gene_data.gtex <- data.frame(gene_id, bp_length, gene_symbol)
gene_data.gtex <- distinct(gene_data.gtex, gene_symbol, .keep_all = TRUE)
gene_data.gtex <- gene_data.gtex[!is.na(gene_data.gtex$gene_symbol), ]
# Subset the expression data
expr.gtex <- expr.gtex[rownames(expr.gtex) %in% gene_data.gtex$gene_id, ]
expr.gtex <- expr.gtex[match(gene_data.gtex$gene_id, rownames(expr.gtex)), ]
rownames(expr.gtex) <- gene_data.gtex$gene_symbol

saveRDS(expr.gtex, paste0(processed_gtexDir,"/expr_gtex.rds"))
saveRDS(sample_metadata.gtex, paste0(processed_gtexDir,"/sample_metadata.rds"))
saveRDS(gene_data.gtex, paste0(processed_gtexDir,"gene_data.rds"))

