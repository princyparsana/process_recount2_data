homeDir=`pwd | sed -e 's/\/shellscripts//g'`
scriptDir=$homeDir/src
datDir=$homeDir/data
logDir=$homeDir/log
rpkm=$homeDir/data/rpkm
auc=$homeDir/data/auc_scaled
raw=$homeDir/data/raw
srcDir=$homeDir/src/

mkdir $datDir
mkdir $logDir
mkdir $rpkm
mkdir $auc
mkdir $raw

cd $datDir
Rscript $srcDir/downloadData.R $raw $auc $rpkm >$logDir/downloadData.log
Rscript $srcDir/processing_smallrna_replicates.R >$logDir/process_smallrna_replicates.log
