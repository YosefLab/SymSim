library("devtools")
load_all("SymSim")
library('plyr')
load_all('easyGgplot2')

args = commandArgs(trailingOnly=TRUE)
ncols <- 2000
colfunc <- colorRampPalette(c("blue", "green", "red", "purple"))
allcol <- colfunc(ncols)

#####################################################################
####### Loading Old Simulation
#####################################################################
dir <- args[1]
outputdir <- args[2]
filenames <- list.files(dir)
lapply(filenames,function(filename){

load(paste(dir,filename,sep=''))
k <- strsplit(filename,'_')[[1]][2]; k <- strsplit(k,'.',fixed=T)[[1]][1];k=as.numeric(k)
sim_params[k,]
protocol <- sim_params$protocol[k]

if(protocol=='10x'){
		load('expression_mRNA_17-Aug-2014.robj')
		counts <- counts[,as.numeric(meta[1,])==3 & !is.na(as.numeric(meta[1,])==3)]
		Cortex_counts <- counts
}else if(protocol=='ss2'){
	load('ExperimentalData/130cells.raw.counts.robj')
	Tcells_counts <- counts
}

true_counts <- true_counts_res[[1]]
obs_counts <- observed_counts[[1]]
zero_rows <- rowSums(obs_counts)==0
zero_cols <- colSums(obs_counts)==0


#####################################################################
####### Plotting heatmap
#####################################################################

p1=PlotCountHeatmap(LogDist(true_counts,seq(0, 4, 0.4)),rowMeans(true_counts),
given_ord= NA,zeropropthres=1,filename=NA,data_name='true_counts',saving=F)

if(protocol=='10x'){plotname <- "UMI sim"}else if(protocol=='ss2'){plotname <-"SS2 sim" }

p2=PlotCountHeatmap(LogDist(obs_counts,seq(0, 4, 0.4)),rowMeans(obs_counts),
                 NA,1, data_name = plotname)

# p4=PlotCountHeatmap(LogDist(Tcells_counts,seq(0, 4, 0.4)),rowMeans(Tcells_counts),
#                  NA,1, data_name = "ss2 T-cell data")
# p5=PlotCountHeatmap(LogDist(Cortex_counts,seq(0, 4, 0.4)),rowMeans(Cortex_counts ),
#                  NA,1, data_name = "UMI Cortex data")

# # save(p4,p5,file='ExperimentalData/realdata_exprs_heatmap.robj')

load('ExperimentalData/realdata_exprs_heatmap.robj')
load("ExperimentalData/realdata_expr_matrices.RData")

outputname <- paste(outputdir,'figure4.param.',k,'.pdf',sep='')
pdf(outputname)
if(protocol=='10x'){
ggplot2.multiplot(p1[[2]], p2[[2]],p5[[2]], cols=2)
}else if(protocol=='ss2'){
ggplot2.multiplot(p1[[2]], p2[[2]],p4[[2]], cols=2)
}
#####################################################################
####### Plotting FNR curves
#####################################################################
par(mfrow=c(3,3))
plotFNRsim(current_counts =obs_counts[!zero_rows,!zero_cols],
 true_counts = true_counts[!zero_rows,!zero_cols],
 x_use = "true", titlestr = paste(plotname,'x true counts'))
plotFNRsim(current_counts = obs_counts[!zero_rows,!zero_cols],
 true_counts = true_counts[!zero_rows,!zero_cols],
  x_use = "observed", titlestr = paste(plotname,'x observed counts'))

if(protocol=='10x'){
	exp_counts <- Cortex_counts
	expname <- '10x Cortex data'
}else if(protocol=='ss2'){
	exp_counts <- expr_matrix_ss2_4heatmap
		expname <- 'SS2 Tcell data'
}

ngenes <- dim(exp_counts)[1]
nhk <- ngenes*0.15
cv_vec <- apply(exp_counts, 1, cv)
mean_vec <- rowMeans(log2(exp_counts+1))
mean_vec[which(cv_vec>0.5)] <- 0
temp <- sort.int(mean_vec, decreasing = T, index.return = T)
hk <- temp$ix[1:nhk]
col_vec_ss2 <- allcol[round(rescale2range(colSums(exp_counts), 2000))]
plotFNRreal(exp_counts, hk, data_name=expname, col_vec_ss2)


#####################################################################
####### Plotting QQ plot
#####################################################################

# plot comparison of distributions of means, 0%, fano


percent_nonzero <- function(x) {return(sum(x>0)/length(x))}

qqout <- qqplot(rowMeans(log2(obs_counts+1)), rowMeans(log2(exp_counts+1)), plot.it = TRUE, 
                col=adjustcolor("blue", alpha.f = 0.5), xlab = "simulated observed counts", ylab = "real data counts",
                main="mean of log2 counts", cex.lab=1.5)
abline(c(0,0),c(1,1),col='red')
qqout <- qqplot(apply(obs_counts, 1, percent_nonzero), apply(exp_counts, 1, percent_nonzero), plot.it = TRUE, 
                col=adjustcolor("blue", alpha.f = 0.5), xlab = "simulated observed data", ylab = "real data",
                main="percent non zero", cex.lab=1.5)
abline(c(0,0),c(1,1),col='red')
qqout <- qqplot(apply(log2(obs_counts+1), 1, fano), apply(log2(exp_counts+1), 1, fano), plot.it = TRUE, 
                col=adjustcolor("blue", alpha.f = 0.5), xlab = "simulated observed data", ylab = "real data",
                main="fano of log2 counts", cex.lab=1.5)
abline(c(0,0),c(1,1),col='red')


#####################################################################
####### Plotting cell size distribution
#####################################################################

hist(colSums(obs_counts),col='coral')
hist(colSums(exp_counts)/length(exp_counts[,1])*10000,col='coral4')

dev.off()
})

