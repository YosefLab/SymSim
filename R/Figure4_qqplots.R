
## umi
counts <- cortex_counts[rowSums(cortex_counts>1)>10, ]
best_matches <- BestMatchParams('umi',counts,'best_params.umi.qqplot.pdf', n_optimal=10)
write.table(best_matches, "best_match_params_umi.txt", row.names = F,  quote=F, sep="\t")

best_k <- as.numeric(rownames(best_matches)[1])
load("SymSim/grid_summary/exp_figure4umi_Lgrid.summary.robj")

mean_exprs <- quantile(rowMeans(counts+1,na.rm=T),seq(0,1,0.002))
fano_exprs <- quantile(apply(counts,1,fano),seq(0,1,0.002),na.rm=T)
percent0 <- quantile(apply(counts,1,percent_nonzero),seq(0,1,0.002))

grid_summary <- list(mean_bins,fano_bins,nonzero_bins)
exp_summary <- list(mean_exprs,fano_exprs,percent0)
umi_data <- lapply(c(1:3),function(k){
  bin1=grid_summary[[k]][best_k,]
  bin2=exp_summary[[k]]
  if(k%in%c(1,3)){
    bin1=log(base=10,bin1);bin2=log(base=10,bin2)
  }
  data.frame(simulation=bin1,experimental=bin2,summary=c('Mean','Fano','Percent Nonzero')[k],tech='umi')
})

## ss2
load("ExperimentalData/Th17_data.RData")
load("ExperimentalData/expressed_genes_ss2.RData")
expressed_genes <- which(rowSums(ss2_130cells_counts > 1) > 10)
exp_counts <- ss2_130cells_counts[expressed_genes,]
counts <- exp_counts

best_matches <- BestMatchParams('ss2',counts,'best_params.ss2.qqplot.pdf', n_optimal=10)
write.table(best_matches, "best_match_params_ss2.txt", row.names = F,  quote=F, sep="\t")

best_k <- as.numeric(rownames(best_matches)[1])
load("SymSim/grid_summary/exp_figure4ss2_Lgrid.summary.robj")

mean_exprs <- quantile(rowMeans(counts+1,na.rm=T),seq(0,1,0.002))
fano_exprs <- quantile(apply(counts,1,fano),seq(0,1,0.002),na.rm=T)
percent0 <- quantile(apply(counts,1,percent_nonzero),seq(0,1,0.002))

grid_summary <- list(mean_bins,fano_bins,nonzero_bins)
exp_summary <- list(mean_exprs,fano_exprs,percent0)
ss2_data <- lapply(c(1:3),function(k){
  bin1=grid_summary[[k]][best_k,]
  bin2=exp_summary[[k]]
  if(k%in%c(1,3)){
    bin1=log(base=10,bin1);bin2=log(base=10,bin2)
  }
  data.frame(simulation=bin1,experimental=bin2,summary=c('Mean','Fano','Percent Nonzero')[k],tech='ss2')
})

plot_data <- c(ss2_data,umi_data)
p=lapply(plot_data,function(X){
  ggplot(X, aes(x=experimental, y=simulation)) +
    geom_point() + geom_abline(intercept = 0, slope = 1,col='red')
})
multiplot(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]], cols=3)


pdf(file=plotfilename)
par(mfrow=c(3,3))
for(i in c(1:3)){
  # for(j in c(1:3)){
  for(k in c(1:3)){
    bin1=grid_summary[[k]][best_match[i],]
    bin2=exp_summary[[k]]
    if(k %in% c(1,3)){bin1 <- log(base=10,bin1);bin2 <- log(base=10,bin2)}
    plot(bin1,bin2,pch=16,xlab='simulated values',ylab='experimental values',main=paste(i,'best','match', plotnames[k]))
    lines(c(-10,10),c(-10,10),col='red')      
  }
}
# }
dev.off()

plot_data <- c(ss2_data,umi_data)
p=lapply(plot_data,function(X){
  ggplot(X, aes(x=experimental, y=simulation)) +
    geom_point() + geom_abline(intercept = 0, slope = 1,col='red')
})
