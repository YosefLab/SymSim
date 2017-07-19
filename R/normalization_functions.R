TestScaling <- function(sim_results){
	true <- sim_results[[2]]
	sampled <- sim_results[[3]]
	biased <- sim_results[[4]]
	norm1 <- Normalize(sampled)
	norm2 <- Normalize(biased)
	cor_compare <- lapply(c(1:length(true[,1])),function(i){
		cor1<-c(cor(true[i,],sampled[i,]),sapply(norm1[[2]],function(X){cor(true[i,],X[i,])}))
		cor2<-c(cor(true[i,],biased[i,]), sapply(norm2[[2]],function(X){cor(true[i,],X[i,])}))
		return(c(cor1,cor2))
	})
	cor_compare <- do.call(rbind,cor_compare)
	return(cor_compare)
}

#' Run all normalization
#' 
#' This function runs rpm, deseq, tmm, uq (quantile normalization), scran, and clt normalization
#' @param x matrix of transcript counts where rows are genes and columns are cells
#' @return list(sf,results) sf is a nmethods * ncells matrix that stores the scaling factors of different methods, results are the gene-cell matrix after normalization
library('DESeq2')
library('scran')
library('scone')
library('edgeR')


Normalize <- function(x){
	temp_t <- proc.time()[3]
	x_rpm <- rpm(x)
	print(paste('rpm',proc.time()[3]-temp_t))
	temp_t <- proc.time()[3]
	x_deseq <- deseq(x)
	print(paste('deseq',proc.time()[3]-temp_t))
	temp_t <- proc.time()[3]
	x_tmm <- tmm(x)
	print(paste('tmm',proc.time()[3]-temp_t))
	temp_t <- proc.time()[3]
	x_uq <- uq(x)
	print(paste('uq',proc.time()[3]-temp_t))
	temp_t <- proc.time()[3]
	if(length(x[1,])>1000){x_scran <- scran_norm(x)}else{x_scran=list(NA,NA)}
	# x_census <- census_norm(x)
	x_clt <- CLT_norm(x)
	print(paste('clt',proc.time()[3]-temp_t))
	temp_t <- proc.time()[3]
	results <- list(x_rpm,x_deseq,x_tmm,x_uq,x_scran,x_clt)
	names(results) <- c('rpm','deseq','tmm','uq','scran','clt')
	sf <- lapply(results,function(X){X[[1]]})
	results <- results[!is.na(sf)]
	sf <- sf[!is.na(sf)]
	sf <- do.call(rbind,sf)
	rownames(sf) <- names(results)
	results <- lapply(results,function(X){X[[2]]})
	return(list(sf,results))
}


#' RPM normalization
#' 
#' This function runs the Reads Per Million normalization (from Vallejos2017 supplementary material), so the total read counts per cell is made equal across all cells
#' @param x matrix of transcript counts where rows are genes and columns are cells
#' @return list(s,counts) s is a ncells long vector that stores the scaling factors for each cell, counts is the gene-cell matrix after normalization
rpm <- function(x) {
	s <- colSums(x)
	n.cells<-ncol(x)
	s <- n.cells*(s/sum(s))
	counts <- t(t(x)/s)
	return(list(s=s, counts=counts))
}

#' RPM normalization
#' 
#' This function runs the DESeq normalization (from Vallejos2017 supplementary material), DESeq defines scaling factor estimates based on a pseudoreference sample, which is built with the geometric mean of gene counts across all cells.
#' @param x matrix of transcript counts where rows are genes and columns are cells
#' @return list(s,counts) s is a ncells long vector that stores the scaling factors for each cell, counts is the gene-cell matrix after normalization
deseq <- function(x) {
	geoMeans <- apply(x, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))) #set the geomMean method we want```
	s <-estimateSizeFactorsForMatrix(x, geoMeans=geoMeans) #estimate size factors```
	n.cells<-ncol(x)
	s <- n.cells*(s/sum(s))
	counts <- t(t(x)/s)
	return(list(s=s, counts=counts))
}

#' TMM normalization
#' 
#' This function runs the Trimmed Mean of M values (TMM) normalization (from Vallejos2017 supplementary material), TMM trims away extreme log fold changes to normalize the counts based on the remaining set of nondifferentially expressed genes.
#' @param x matrix of transcript counts where rows are genes and columns are cells
#' @return list(s,counts) s is a ncells long vector that stores the scaling factors for each cell, counts is the gene-cell matrix after normalization
tmm <- function(x) {
	s <- calcNormFactors(x)*colSums(x)
	n.cells<-ncol(x)
	s <- n.cells*(s/sum(s))
	counts <- t(t(x)/s)
	return(list(s=s, counts=counts))
}
#' UQ normalization
#' 
#' This function runs the Upper Quantile (UQ) normalization (from Vallejos2017 supplementary material), scaling factor estimates is proportional to the 75th percentile of the distribution of counts within each cell
#' @param x matrix of transcript counts where rows are genes and columns are cells
#' @return list(s,counts) s is a ncells long vector that stores the scaling factors for each cell, counts is the gene-cell matrix after normalization
uq <- function(x) {
	s <- apply(x, 2, quantile, .75)
	n.cells<-ncol(x)
	s <- n.cells*(s/sum(s))
	counts <- t(t(x)/s)
	return(list(s=s, counts=counts))
}
