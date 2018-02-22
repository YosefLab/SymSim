library("devtools")
library('plyr')
load_all("SymSim")
load_all('easyGgplot2')


load("Rvariables/sim_params_1popFigure4_0212.RData")

percent_nonzero <- function(x) {return(sum(x>0)/length(x))}

protocol <- "umi"
load('ExperimentalData/expression_mRNA_17-Aug-2014.robj') # load cortex data before imputation
cortex_counts <- counts[,as.numeric(meta[1,])==3 & !is.na(as.numeric(meta[1,])==3)]
rm(counts)
p1=PlotCountHeatmap(LogDist(cortex_counts,seq(0, 4, 0.4)),rowMeans(cortex_counts),
                    given_ord= NA,zeropropthres=1,filename=NA,data_name='cortex_counts_obs',saving=F)
data_train <- read.table("ExperimentalData/data_train",  header = F, stringsAsFactors = F)
label_train <- read.table("ExperimentalData/label_train", header=F, stringsAsFactors = F)
cortex_counts_imputed <- t(as.matrix(data_train[label_train$V1==3,]))
p2=PlotCountHeatmap(LogDist(cortex_counts_imputed,seq(0, 4, 0.4)),rowMeans(cortex_counts_imputed),
                    given_ord= NA,zeropropthres=1,filename=NA,data_name='cortex_counts_imputed',saving=F)


exp_counts <- cortex_counts[rowSums(cortex_counts)>0, ]
#exp_counts <- t(t(exp_counts)/colSums(exp_counts)) * 10^6
nparams <- dim(sim_params)[1]
mean_dist <- numeric(nparams); nonzero_dist <- numeric(nparams); fano_dist <- numeric(nparams)
for (k in 213:dim(sim_params)[1]){
  datafn <- sprintf("Rvariables/exp_figure4umi10/output_%d_true.RData",k); load(datafn)
  datafn <- sprintf("Rvariables/exp_figure4umi10/output_%d_obs.RData",k); load(datafn)
  
  true_counts <- true_counts_res[[1]]
  #true_counts <- t(t(true_counts)/colSums(true_counts)) * 10^6
  obs_counts <- observed_counts[[1]][rowSums(observed_counts[[1]])>0, ]
  #obs_counts <- t(t(obs_counts)/colSums(obs_counts)) * 10^6
  
  pdf(sprintf("plots/figure4_0214/umi10_%d_alpha%4.3f_geffect_sd%3.1f_depth%d_Sigma%2.1f.pdf", k, 
              sim_params$alpha_mean[k], sim_params$gene_effects_sd[k], sim_params$depth_mean[k],
              sim_params$Sigma[k]), 9, 12)
  
  # p3=PlotCountHeatmap(LogDist(true_counts,seq(0, 4, 0.4)),rowMeans(true_counts),
  #                     given_ord= NA,zeropropthres=1,filename=NA,data_name='true_counts',saving=F)
  # if(protocol=='umi'){plotname <- "UMI sim"}else if(protocol=='ss2'){plotname <-"SS2 sim" }
  # p4=PlotCountHeatmap(LogDist(obs_counts,seq(0, 4, 0.4)),rowMeans(obs_counts),
  #                     NA,1, data_name = plotname)
  # ggplot2.multiplot(p1[[2]], p2[[2]], p3[[2]], p4[[2]], cols=2)
  
  p_nonzero_true <- apply(true_counts, 1, percent_nonzero)
  p_nonzero_obs <- apply(obs_counts, 1, percent_nonzero)
  p_nonzero_exp <- apply(exp_counts, 1, percent_nonzero)
  p_nonzero_exp_imputed <- apply(cortex_counts_imputed, 1, percent_nonzero)
  mean_true <- rowMeans(log2(true_counts+1))
  mean_obs <- rowMeans(log2(obs_counts+1))
  mean_exp <- rowMeans(log2(exp_counts+1))
  fano_obs <- apply(log2(obs_counts+1), 1, fano)
  fano_exp <- apply(log2(exp_counts+1), 1, fano)
  
  minlen <- min(length(obs_counts[,1]),length(exp_counts[,1]))
  mean_dist[k] <- ks.dist(sort(sample(mean_exp, minlen, replace = F)), sort(sample(mean_obs, minlen, replace = F)))$D
  nonzero_dist[k] <- ks.dist(sort(sample(p_nonzero_obs, minlen, replace = F)), sort(sample(p_nonzero_exp, minlen, replace = F)))$D
  fano_dist[k] <- ks.dist(sort(sample(fano_obs, minlen, replace = F)), sort(sample(fano_exp, minlen, replace = F)))$D
  
  par(mfrow=c(4,3), cex.lab=1.2)
  hist(p_nonzero_true, 50, main="simulated true", ylim = c(0, 7000), xlab="percent non-zero")
  hist(p_nonzero_obs, 50, main="simulated obs", ylim = c(0, 7000), xlab="percent non-zero", col = "green")
  hist(p_nonzero_exp, 50, main="real data", ylim = c(0, 7000), xlab="percent non-zero", col="red")
  
  hist(mean_true,50, main="simulated true", ylim = c(0, 6000), xlab="mean of log2 counts")
  hist(mean_obs, seq(0,max(mean_exp)+0.5,0.2), xlim=c(0,max(mean_exp)), main="simulated obs", ylim = c(0, 6000), xlab="mean of log2 counts", col = "green")
  hist(mean_exp, seq(0,max(mean_exp)+0.5,0.2), xlim=c(0,max(mean_exp)), main="real data", ylim = c(0, 6000), xlab="mean of log2 counts", col="red")
  
  qqout <- qqplot(mean_obs, mean_exp, plot.it = TRUE, col=adjustcolor("blue", alpha.f = 0.5), 
                  xlab = "simulated observed counts", ylab = "real data counts", 
                  main=sprintf("mean(log2(counts)), ks.dist=%4.3f", mean_dist[k]))
  abline(c(0,0),c(1,1),col='red')
  
  qqout <- qqplot(p_nonzero_obs, p_nonzero_exp, plot.it = TRUE, col=adjustcolor("blue", alpha.f = 0.5), 
                  xlab = "simulated observed data", ylab = "real data", 
                  main=sprintf("percent >0, ks.dist=%4.3f", nonzero_dist[k]))
  abline(c(0,0),c(1,1),col='red')
  
  qqout <- qqplot(p_nonzero_true, p_nonzero_exp_imputed, plot.it = TRUE, col=adjustcolor("blue", alpha.f = 0.5), 
                  xlab = "simulated observed data", ylab = "real data imputed", main="percent non zero")
  abline(c(0,0),c(1,1),col='red')
  
  qqout <- qqplot(fano_obs, fano_exp, plot.it = TRUE,  col=adjustcolor("blue", alpha.f = 0.5), 
                  xlab = "simulated observed data", ylab = "real data", 
                  main=sprintf("fano(log2(counts)), ks.dist=%4.3f", fano_dist[k]))
  abline(c(0,0),c(1,1),col='red')
  
  hist(colSums(obs_counts),col='coral', main=sprintf("obs cell size %7.1f", mean(colSums(obs_counts))))
  hist(colSums(exp_counts),col='coral4', main=sprintf("exp cell size %7.1f", mean(colSums(exp_counts))))
  
  dev.off()
}




qqplot(p_nonzero_obs, p_nonzero_obs*2, plot.it = T)

# pdf("test.pdf")
# p1[[2]]
# p2[[2]]
# par(mfrow=c(1,2))
# plot(1,1)
# plot(2,2)
# dev.off()
###############################
load_all("SymSim")
true_counts_res <- SimulateTrueCounts(ncells_total=50, ngenes=100,
                                      evf_center=1,nevf=10,evf_type="one.population",Sigma=sim_params$Sigma[k], 
                                      gene_effects_sd=sim_params$gene_effects_sd[k],
                                      gene_effect_prob=0.3,bimod=1,randseed=0, param_realdata = "zeisel.imputed", SE=F)

load("SymSim/match_params.zeisel_imputed.robj")
dim(match_params)


load("Rvariables/sim_params_multipopDE.RData")
head(sim_paramsDE)

load("Rvariables/sim_params_1popFigure4_umi_lowalpha.RData")
sim_params$protocol[which(sim_params$protocol=="10x")] <- "umi"
save(list=c("sim_params"), file="Rvariables/sim_params_1popFigure4_umi_lowalpha.RData")

icell <- 1
ngenes <- 20000
ncells <- 1000
observed_counts <- matrix(0, ngenes, ncells)
test <- amplify_1cell(true_counts_1cell =  true_counts_res[[1]][, icell], protocol="umi", 
                      rate_2cap=0.01, gene_len=gene_len, 
                      amp_bias = cal_amp_bias(0.01, nbins=20, gene_len=gene_len, amp_bias_limit=c(-0.2, 0.2)),  
                      rate_2PCR=0.7, nPCR=18, N_molecules_SEQ = 70000) 

hist(apply(exp_counts, 2, cv))
hist(apply(obs_counts, 2, cv))
