#############input arguments###############################
#inputargs <- commandArgs(TRUE)
# dirName <- inputargs[1]
#rawExp <- inputargs[2]
# saveDir <- inputargs[3]
dirName <- "/work-zfs/abattle4/parsana/process_recount2_data/data/sra_metadata/"
rawExp <- "/work-zfs/abattle4/parsana/process_recount2_data/data/rpkm/sra.Rds"
saveDir <- "/work-zfs/abattle4/parsana/process_recount2_data/data/replicates_merged_sra/"
dir.create(saveDir)

library(SRAdb)
library(dplyr)
# load(rawExp)# load data
rse_gene <- readRDS(rawExp)
runids <- rse_gene$run # get run ids

# get corresponding metadata from SRAdb
dbFilename <- paste(dirName, "SRAmetadb.sqlite", sep = "")
#sqlfile <- getSRAdbFile(, destfile = paste(dbFilename, ".gz", sep = "")) ## not working on marcc right now, using the previously downloaded version of db
sra_con <- dbConnect(SQLite(),dbFilename)
res <- dbGetQuery(sra_con, paste("select * from sra where run_accession in ('",paste(runids, collapse = "','"),"')", sep = ""))

run_diff <- setdiff(runids, res$run_accession)
study_diff <- table(rse_gene$project[rse_gene$run %in% run_diff])

## only select runs that have annotations/metadata
rse_gene <- rse_gene[,match(res$run_accession,rse_gene$run)]

potential_cell_lines <- c(grep("cell_line", res$sample_attribute), grep("cells", res$experiment_title))
# res$experiment_title[potential_cell_lines]

smallrna_samples <- res[res$library_selection %in% "size fractionation",]$run_accession
rse_gene <- rse_gene[,!rse_gene$run %in% smallrna_samples]

studies <- unique(rse_gene$project)
# replicate_exp <- table(rse_gene$experiment)>1
# studies <- unique(rse_gene$project[rse_gene$experiment %in% names(replicate_exp)[which(replicate_exp)]])

each_study_count <- lapply(studies, function(x,y){
  print(x)
  this_y <- y[, y$project %in% x]
  this_study <- as.data.frame(t(this_y@assays$data$counts))
  this_study$exp <- this_y$experiment
  if(length(unique(this_y$experiment)) == length(unique(this_y$run))){
    this_count <- this_study
  }else{
    this_count <- this_study %>% group_by(exp) %>% summarise_all(median) %>% as.data.frame
  }
  this_count
}, rse_gene)

each_study_count <- lapply(each_study_count, function(this_study){
  rownames(this_study) <- this_study$exp
  this_study$exp <- NULL
  this_study
})

names(each_study_count) <- studies

## Exclude samples with more than 50% of genes with 0 expression value
total_genes = ncol(each_study_count[[1]])

prop_genes_notexp = lapply(each_study_count, function(x) apply(x, 1, function(y) (sum(y<=0.5)/total_genes)<0.5))

each_study_count = mapply(function(x,y){
  x[y,]
}, each_study_count, prop_genes_notexp, SIMPLIFY = FALSE)

saveRDS(each_study_count, file = paste(saveDir,"rpkm_replicates_merged.Rds", sep = ""))
