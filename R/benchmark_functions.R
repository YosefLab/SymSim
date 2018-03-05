#' Perform tSNE 
#' 
#' This function takes the ground truth (meta), the expresssion matrix (data), and performs tSNE, saves the result to a jpeg file with user specified name
#' @param meta: simulation parameters
#' @param data: expression matrix
#' @param plotname: the name of the jpeg file
#' @param label: the column name of the meta data that the points needs to be colored by

<<<<<<< HEAD
PlotTsne <- function(meta,data,plotname,label,n_pc=20,evf_type="discrete",saving=F,perplexity=10){
=======
PlotTsne <- function(meta,data,plotname,label,evf_type="discrete",saving=F,perplexity=10){
>>>>>>> b1c727f6ac83412053a415eccd221e3a5179bbd7
  library('Rtsne')
  uniqcols<-c(1:length(data[1,]))[!duplicated(t(data))]
  data <- data[,uniqcols];meta <- meta[uniqcols,,drop=FALSE] 
  uniqrows<-c(1:length(data[,1]))[!duplicated(data)]
  data <- data[uniqrows,]
<<<<<<< HEAD
  data_pc <- prcomp(t(data))
  data_pc <- data_pc$x[,c(1:n_pc)]
  data_tsne=Rtsne(data_pc,perplexity=perplexity)
=======
  data_tsne=Rtsne(t(data),perplexity=perplexity)
>>>>>>> b1c727f6ac83412053a415eccd221e3a5179bbd7
  
  if(evf_type=="discrete" | evf_type=="one.population"){
    plot_tsne <- cbind(meta, data.frame(label=factor(meta[,label]),x=data_tsne$Y[,1],y=data_tsne$Y[,2]))
  }else{
    plot_tsne <- cbind(meta, data.frame(label=meta[,label],x=data_tsne$Y[,1],y=data_tsne$Y[,2]))
  }
  p <- ggplot(plot_tsne, aes(x, y))
  p <- p + geom_point()
  p <- p + geom_point(aes(colour = plot_tsne[['label']]))+labs(color=label)
  if(saving==T){ggsave(p,filename=plotname,device='pdf',width=5,height=4)}
  if(saving==F){p <- p + ggtitle(plotname)}
  return(list(plot_tsne,p))
}

#' Perform PCA  
#' 
#' This function takes the ground truth (meta), the expresssion matrix (data), and performs tSNE, saves the result to a jpeg file with user specified name
#' @param meta: simulation parameters
#' @param data: expression matrix
#' @param plotname: the name of the jpeg file
#' @param label: the column name of the meta data that the points needs to be colored by

PlotPCA <- function(meta,data,plotname,label,discrete=T,saving=F){
    uniqcols<-c(1:length(data[1,]))[!duplicated(t(data))]
    data <- data[,uniqcols];meta <- meta[uniqcols,] 
    uniqrows<-c(1:length(data[,1]))[!duplicated(data)]
    data <- data[uniqrows,]
	data_pc=prcomp(t(data))
	if(discrete==T){
		plot_pca=data.frame(meta,label=factor(meta[,label]),x=data_pc$x[,1],y=data_pc$x[,2])		
	}else{
		plot_pca=data.frame(meta,label=meta[,label],x=data_pc$x[,1],y=data_pc$x[,2])
	}
	p <- ggplot(plot_pca, aes(x, y))
	p <- p + geom_point()
	p <- p + geom_point(aes(colour = plot_pca[['label']]))+labs(color=label)
	if(saving==T){ggsave(p,filename=plotname,device='jpeg',width=5,height=4)}
	return(list(plot_pca,p))
}

DEseq_test <- function(result){
	evfs <- result[[2]][,grep('evf',colnames(result[[2]]))]
	pop1 <- result[[2]]$pop==1
	de_true <- lapply(c('kon_effect','koff_effect','s_effect'),function(X){
	    effect <- result[[3]][,grep(X,colnames(result[[3]]))]
	    value <- t(as.matrix(evfs) %*% t(as.matrix(effect)))
	    diff <- lapply(c(1:ngenes),function(gid){
	        diff <- t.test(value[gid,pop1],value[gid,!pop1])
	        c(diff$p.value,diff$conf.int)
	    })
	    diff <- do.call(rbind,diff)
	    pos <- c(1:ngenes)[diff[,1]<0.05 & diff[,2] > 0.1]
	    neg <- c(1:ngenes)[diff[,1]<0.05 & diff[,3] < (-0.1)]
	    pos <- pos[!is.na(pos)];neg <- neg[!is.na(neg)]
        return(list(pos,neg))
	})
	logFC <- log2(rowMeans(result[[1]][[1]][,pop1])/ rowMeans(result[[1]][[1]][,!pop1]))
	categories <- c('kon_pos','kon_neg','koff_pos','koff_neg','s_pos','s_neg')
	colnames(result[[1]][[3]]) <- paste0("cell", 1:length(result[[1]][[3]][1,]))
	rownames(result[[1]][[3]]) <- paste0("gene", 1:length(result[[1]][[3]][,1]))
	coldata <- data.frame(pop=paste0("pop", result[[2]]$pop==1))
	rownames(coldata) <- colnames(result[[1]][[3]])
	dds <- DESeqDataSetFromMatrix(countData = round(result[[1]][[3]])+1,colData = coldata, design= ~ pop)
	dds <- DESeq(dds)
	rnms <- resultsNames(dds)
	res <- results(dds)
	DEseq_pos <- c(1:ngenes)[res$pvalue<0.05 & res$log2FoldChange>0]
	DEseq_neg <- c(1:ngenes)[res$pvalue<0.05 & res$log2FoldChange<0]
	overlaps <- do.call(rbind,
	    lapply(do.call(c,de_true),function(X){
	        c(length(intersect(X,DEseq_neg))/length(X),length(intersect(X,DEseq_pos))/length(X))
	    })
	)
	rownames(overlaps) <- categories
	colnames(overlaps) <- c('DE_pos','DE_neg')
	return(list(logFC,res$log2FoldChange,res$pvalue,overlaps))
}