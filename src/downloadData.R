select.genes <- function(rse.object, threshold, ...){
  counts <- SummarizedExperiment::assay(rse.object, 1)
  min.samples <- round(0.25*dim(rse.object)[2]) # default min.samples expression in at least 1/4 of samples
  keep <- apply(counts, 1, function(x, n = min.samples){
    t = sum(x >= threshold) >= n
    t
  })
  rse.object <- rse.object[keep,]
}

###########################################

pc.genes <- read.delim("/work-zfs/abattle4/parsana/networks_correction/data/etc/protein_coding.txt", 
header = F, stringsAsFactors = F)
pc.genes <- pc.genes$V2[!pc.genes$V1 %in% c("chrM","chrY")]

overlapping_genes <- read.delim("/work-zfs/abattle4/parsana/networks_correction/data/etc/ensembl_ids_overlapping_genes.txt",
        stringsAsFactors = F)


studies_interest <- c("sra", "TCGA", "SRP012682")
sapply(studies_interest, download_study)
sapply(studies_interest, function(eachstudy, pc_genes, ov_genes){
        load(file.path(eachstudy, "rse_gene.Rdata"))
        rse_gene <- scale_counts(rse_gene, by = "auc", round = FALSE)
        rse_gene <- select.genes(rse_gene, threshold = 0.1)
        rse_gene <- rse_gene[rownames(rse_gene) %in% pc_genes,]
        rse_gene <- rse_gene[!rownames(rse_gene) %in% ov_genes,]
        counts <- SummarizedExperiment::assay(rse_gene, 1)
        SummarizedExperiment::assay(rse_gene, 1) <- log2(counts+2)
        saveRDS(rse_gene, file = paste("/work-zfs/abattle4/parsana/recount_networks/data/", eachstudy,".Rds", sep = ""))
        }, pc.genes, overlapping_genes$x)


