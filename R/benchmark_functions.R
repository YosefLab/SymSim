#' cluster-ability of a dataset 
#' 
#'how well can the cell population be reconstructed by dimensionality reduction methods
#' @param data: transcriptomics matrix
#' @param meta: information about cells, must contain 'pop' column
#' @param n_pc: number of principal components that is fed to tsne
#' @param perplexity: perplexity parameter for tsne
#' @param dims: number of dimensions kept after tsne

ClusterQuality <- function(data,meta,n_pc,perplexity,dims,npop=5,return_kmeans=F,pca_scale=F){
  uniqcols<-c(1:length(data[1,]))[!duplicated(t(data))]
  data <- data[,uniqcols];meta <- meta[uniqcols,,drop=FALSE] 
  uniqrows<-c(1:length(data[,1]))[!duplicated(data)]
  data <- data[uniqrows,]
  data_pc <- prcomp(t(data), scale. = pca_scale)
  data_pc <- data_pc$x[,c(1:n_pc)]
  data_tsne=Rtsne(dims = dims ,data_pc,perplexity=perplexity)
  lowdim_data <- data_tsne$Y
  assignment <- kmeans(lowdim_data,npop,nstart=20)
  ri <- rand.index(assignment[[1]],meta$pop)
  lowdim_dist <- as.matrix(dist(lowdim_data, method = "euclidean"))
  sw <- silhouette(x=meta$pop, dist=dist(lowdim_data, method = "euclidean"))
  sw_summ <- summary(sw)
  sw_mean <- sw_summ[['avg.width']]
  minpop <- sw_summ[["clus.sizes"]]==min(sw_summ[["clus.sizes"]])
  sw_minpop <- sw_summ[["clus.avg.widths"]][minpop]
  # for all cells that are in pop2, what is the most common cluster they are in
  # for all clusters, which one has the highest prop of pop2 cells
  # for the most common pop2 cluster, what is the proportion of pop2 cells? 
  min_in_clst <- sapply(split(x=meta$pop,f=assignment[[1]]),function(X){
    sum(X%in%c(1:npop)[minpop])})
  clst_w_min <- sapply(split(x=meta$pop,f=assignment[[1]]),function(X){
    sum(X%in%c(1:npop)[minpop])/length(X)})
  res <- c(ri,sw_mean,sw_minpop,min_in_clst,clst_w_min)
  names(res) <- c("RI", "SW_mean", "SW_minpop", paste0("min_in_clst",1:npop), paste0("clst_w_min", 1:npop))
  if (return_kmeans)
    return(list(measures=res, kmeans_res=assignment))
  else
    return(res)
}
#' rand index
#'
#' compare two clustering result: proportion of pairs of individual that share the same clustering result in both groupings
#' @param group1 first clustering
#' @param group2 second clustering 

rand.index<-function (group1, group2) 
{
  x <- c(sapply(group1, function(x) {x==group1}))
  y <- c(sapply(group2, function(x) {x==group2}))
  same <- sum(x==y)
  ri <- same/length(x)
  return(ri)
}

#' Perform tSNE 
#' 
#' This function takes the ground truth (meta), the expresssion matrix (data), and performs tSNE, saves the result to a jpeg file with user specified name
#' @param meta: simulation parameters
#' @param data: expression matrix
#' @param plotname: the name of the jpeg file
#' @param label: the column name of the meta data that the points needs to be colored by


PlotTsne <- function(meta, data, evf_type, pca = T, n_pc, perplexity=30, label, saving=F, plotname,system.color=T){
  library('Rtsne')
  uniqcols<-c(1:length(data[1,]))[!duplicated(t(data))]
  data <- data[,uniqcols];meta <- meta[uniqcols,,drop=FALSE] 
  uniqrows<-c(1:length(data[,1]))[!duplicated(data)]
  data <- data[uniqrows,]
  if (pca){
    data_pc <- prcomp(t(data))
    data <- t(data_pc$x[,c(1:n_pc)])
  }
  data_tsne=Rtsne(t(data),perplexity=perplexity)
  plot_tsne <- cbind(meta, data.frame(label=factor(meta[,label]),x=data_tsne$Y[,1],y=data_tsne$Y[,2]))
  p <- ggplot(plot_tsne, aes(x, y))
  p <- p + geom_point(aes(colour = plot_tsne[['label']]),shape=20) + labs(color=label)
  if (system.color==F){
    if(evf_type=="discrete" | evf_type=="one.population"){
      color_5pop <- c("#F67670", "#0101FF", "#005826", "#A646F5", "#980000")
      names(color_5pop) <- 1:5
    }else{
      color_5pop <- c("#CC9521", "#1EBDC8", "#0101FF", "#005826", "#7CCC1F", "#A646F5", "#980000", "#F67670")
      names(color_5pop) <- c("6_7","7_8","8_2","8_3","7_9","9_4","9_5","6_1")
    }
    p <- p + scale_color_manual(values=color_5pop[levels(plot_tsne[['label']])])
  }
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

# calculate area under curve
cal_AUC <- function(x_vec, y_vec){
  xtemp <- x_vec[2:length(x_vec)]-x_vec[1:(length(x_vec)-1)]
  ytemp <- y_vec[2:length(y_vec)]-y_vec[1:(length(y_vec)-1)]
  return(sum(xtemp*(y_vec[1:(length(y_vec)-1)] + ytemp/2)))
}

# calculate sensitivity and specificity
sens_and_spec <- function(isDE_gold, isDE_pred){
  TP <- sum(isDE_gold & isDE_pred, na.rm = T)
  TN <- sum(!isDE_gold & !isDE_pred, na.rm = T)
  FP <- sum(!isDE_gold & isDE_pred, na.rm = T)
  FN <- sum(isDE_gold & !isDE_pred, na.rm = T)
  sensi <- TP/(TP+FN); speci <- TN/(TN+FP)
  return(c(sensi, speci))
}

