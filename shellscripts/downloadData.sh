#!/bin/sh
#SBATCH --time=1:0:0
#SBATCH --mem=40G
#SBATCH --partition=shared


cd /work-zfs/abattle4/parsana/recount_networks/data/
Rscript /work-zfs/abattle4/parsana/recount_networks/src/downloadData.R >/work-zfs/abattle4/parsana/recount_networks/log/downloadData.log
