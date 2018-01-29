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
outputname <- args[2]
filenames <- list.files(dir)
ks_distances <- lapply(filenames,function(filename){
	load(paste(dir,filename,sep=''))
	k <- strsplit(filename,'_')[[1]][2]; k <- strsplit(k,'.',fixed=T)[[1]][1];k=as.numeric(k)
	sim_params[k,]
	protocol <- sim_params$protocol[k]

	if(protocol=='10x'){
		counts <- read.csv('ExperimentalData/data_train',sep=' ',as.is=T);counts <- t(counts)
		label <- read.csv('ExperimentalData/label_train',sep=' ',as.is=T);label <- label[,1]
		counts <- counts[,label==5]
		Cortex_counts <- counts
	}else if(protocol=='ss2'){
		load('ExperimentalData/realdata_exprs_heatmap.robj')
		load("ExperimentalData/realdata_expr_matrices.RData")
		Tcells_counts <- counts
	}

	true_counts <- true_counts_res[[1]]
	obs_counts <- observed_counts[[1]]
	zero_rows <- rowSums(obs_counts)==0
	zero_cols <- colSums(obs_counts)==0

	if(protocol=='10x'){
		exp_counts <- Cortex_counts

		expname <- '10x Cortex data'
	}else if(protocol=='ss2'){
		exp_counts <- expr_matrix_ss2_4heatmap
			expname <- 'SS2 Tcell data'
	}

	# percent_nonzero <- function(x) {return(sum(x>0)/length(x))}
	minlen <- min(length(obs_counts[,1]),length(exp_counts[,1]))
	mean_dist <- ks.dist(rowMeans(log2(obs_counts+1))[c(1:minlen)], rowMeans(log2(exp_counts+1))[c(1:minlen)])$D
	nonzero_dist <- ks.dist(apply(obs_counts, 1, percent_nonzero)[c(1:minlen)], apply(exp_counts, 1, percent_nonzero)[c(1:minlen)])$D
	fano1 <- apply(log2(obs_counts+1), 1, fano)
	fano1 <- fano1[!is.na(fano1)]
	fano2 <-  apply(exp_counts, 1, fano)
	fano2 <- fano2[!is.na(fano2)]
	minlen <- min(length(fano1),length(fano2))-1
	fano_dist <- ks.dist(spec1=fano1[c(1:minlen)],spec2=fano2[c(1:minlen)])$D
	return(c(k,mean_dist,nonzero_dist,fano_dist))
})
ks_distances <- do.call(rbind,ks_distances)
# k_temp <- sapply(filenames,function(filename){
# 	k <- strsplit(filename,'_')[[1]][2]
# 	k <- strsplit(k,'.',fixed=T)[[1]][1]
# 	as.numeric(k)
# })
# ks_distances <- cbind(k_temp,ks_distances)
save(sim_params,ks_distances,file=paste(outputname,'.ks_dist.robj',sep=''))