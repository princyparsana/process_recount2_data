#!/bin/bash -l
#SBATCH
#SBATCH --time=24:0:0
#SBATCH --partition=lrgmem
#SBATCH --mem-per-cpu=200G
#SBATCH --mail-type=end
#SBATCH --mail-user=pravich2@jhu.edu

module load R/3.6.1
module load gcc/5.5.0

cd /work-zfs/abattle4/prashanthi/princy_preprocessing/process_recount2_data/shellscripts/
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src/
datDir=$homeDir/data/
logDir=$homeDir/log
rpkm=$homeDir/data/rpkm/
auc=$homeDir/data/auc_scaled/
raw=$homeDir/data/raw/
srcDir=$homeDir/src/
srametDir=$datDir/sra_metadata/

mkdir $datDir
mkdir $logDir
mkdir $rpkm
mkdir $auc
mkdir $raw
mkdir $srametDir
mkdir $datDir/expr_data

cd $datDir
Rscript $srcDir/downloadData.R $raw $auc $rpkm >$logDir/downloadData.log
Rscript $srcDir/processing_smallrna_replicates.R $srametDir $rpkm/sra.Rds $datDir $rpkm >$logDir/process_smallrna_replicates.log
mkdir $datDir/automated_process_output/
cp $datDir/replicates_merged_sra/rpkm_replicates_merged.Rds $datDir/automated_process_output/sra.Rds
cp $rpkm/SRP012682.Rds $datDir/automated_process_output/gtex.Rds
cp $rpkm/TCGA.Rds $datDir/automated_process_output/tcga.Rds
Rscript $srcDir/process_tcga.R $datDir/automated_process_output $datDir/expr_data
