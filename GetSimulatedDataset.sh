#PBS  -N Npop1Batch
#PBS -q yosef3
#PBS -l mem=2g

cd /data/yosef2/users/chenling/scRNASeq-simulation/
Rscript SCsimR/GetSimulatedDataset.multibatches.R 
