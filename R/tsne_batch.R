library("devtools")
load_all("SymSim"); 
#source("functions.R")
library('plyr')


# load("Symsim/Rvariables/sim_params_1popUMI_Lgrid.RData")

phyla <- Phyla5()
ncells=5000
ngenes=1000
load("SymSim/gene_len_pool.RData")
gene_len <- sample(gene_len_pool[which(gene_len_pool>100)], ngenes, replace = FALSE)

true_count <- SimulateTrueCounts(ncells_total=ncells,min_popsize=50,
  ngenes=ngenes, nevf=10,Sigma=0.1,evf_center=0.1,
  phyla=phyla,param_realdata="zeisel.imputed",
  bimod=0,gene_effects_sd=1,gene_effect_prob=0.3,
  evf_type='discrete',n_de_evf=5,randseed=0,SE=F,
  vary='s',joint=F)

observed_counts <- True2ObservedCounts(
  true_counts=true_count[[1]],meta_cell=true_count[[3]],
  protocol="umi",alpha_mean=0.01,
  lenslope=0.01,nbins=20,gene_len=gene_len,
  amp_bias_limit=c(-0.2, 0.2),rate_2PCR=0.8,nPCR=18,
  depth_mean=4e+04, depth_sd=24000)

results <- BatchTrue2ObservedCounts(true_count,'umi',alpha_mean=0.01,lenslope=0.01,
  batch_effect=0.3,depth_mean=4e+04,depth_sd=24000)

p1 = PlotTsne(meta=true_count[[3]],data=true_count[[1]],plotname='temp',label='pop')
p2 = PlotTsne(meta=observed_counts[[2]],data=observed_counts[[1]],plotname='temp',label='pop')
p3 = PlotTsne(meta=results[[2]],data=results[[1]],plotname='temp',label='pop')
p4 = PlotTsne(meta=results[[2]],data=results[[1]],plotname='temp',label='batch_id')
p5 = PlotTsne(meta=results[[2]][results[[2]][,'batch_id']==1,],data=results[[1]][,results[[2]][,'batch_id']==1],plotname='temp',label='pop')

