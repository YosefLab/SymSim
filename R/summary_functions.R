cv <- function(x) {return(sd(x)/mean(x))}
fano <- function(x) {return((sd(x))^2/mean(x))}
percent_nonzero <- function(x) {return(sum(x>0)/length(x))}

#' Get the best matched parameter 
#'
#' This function matches a real dataset to a database of summary information of simulated datasets, plots a qqplot for user to visualize similarity between their dataset and the simulated dataset, and suggests parameters to use in the simulation
#' @param tech 'nonUMI','UMI' the match database are constructed based on reasonable values for each technology
#' @param counts expression matrix
#' @param plotfilename output name for qqplot
#' @param n_optimal number of top parameter configurations to return
#' @return three set of best matching parameters that was used to simulate the best matching dataset to the experimental dataset
BestMatchParams <- function(tech,counts,plotfilename,n_optimal=3){
  counts <- counts[rowSums(counts>0)>10, ]
  mean_exprs <- quantile(rowMeans(counts+1,na.rm=T),seq(0,1,0.002))
  sd_exprs <- quantile(apply(counts,1,sd),seq(0,1,0.002),na.rm=T)
  percent0 <- quantile(apply(counts,1,percent_nonzero),seq(0,1,0.002))
  if(tech=='nonUMI'){
    load('SymSim/grid_summary/figure4_nonUMI_Lgrid.summary.RData')
  }else if(tech =='UMI'){
    load('SymSim/grid_summary/figure4_UMI_Lgrid.summary.RData')   
  }
  grid_summary <- list(mean_bins,nonzero_bins,sd_bins)

  exp_summary <- list(mean_exprs,percent0,sd_exprs)
  dists <- lapply(c(1:3),function(i){
    if (i %in% c(1,3)){
      dist <- apply(grid_summary[[i]],1,function(X){sum(abs(log10(X)-log10(exp_summary[[i]])))})
    } else{
      dist <- apply(grid_summary[[i]],1,function(X){sum(abs(X-exp_summary[[i]]))})
   }
   return(dist)
  })
  dists <- do.call(cbind,dists)
  dists <- apply(dists,2,function(X){X/mean(X)})
  dists <- rowSums(dists)
  best_match <- sort.int(dists, index.return = T)$ix[1:n_optimal]

  best_params <- lapply(best_match,function(X){sim_params[X,]})
  plotnames <- c('mean','percent_nonzero','sd')
  # pdf(file=plotfilename, 10, 23)
  # par(mfrow=c(6,3))
  # for(i in c(1:6)){
  #   for(k in c(1:3)){
  #     bin1=grid_summary[[k]][best_match[i],]
  #     bin2=exp_summary[[k]]
  #     if(k %in% c(1,3)){bin1 <- log(base=10,bin1);bin2 <- log(base=10,bin2)}
  #     plot(bin1,bin2,pch=16,xlab='simulated values',ylab='experimental values',main=paste(i,'best','match', plotnames[k]))
  #     lines(c(-10,10),c(-10,10),col='red')      
  #   }
  # }
  # dev.off()
  par(mfrow=c(1,3))
  for(k in c(1:3)){
    bin1=grid_summary[[k]][best_match[1],]
    bin2=exp_summary[[k]]
    if(k %in% c(1,3)){bin1 <- log(base=10,bin1);bin2 <- log(base=10,bin2)}
    plot(bin1,bin2,pch=16,xlab='simulated values',ylab='experimental values',main=plotnames[k])
    lines(c(-10,10),c(-10,10),col='red')      
  }
  best_params <- do.call(rbind,best_params)
  best_params <- best_params[, c("gene_effects_sd", "gene_effect_prob", "nevf", "Sigma",
                                 "alpha_mean", "alpha_sd", "depth_mean", "depth_sd")]
  return(best_params)
}

#' Get the logged distribution from master equation simulations
#'
#' This function converts the frequency on integers from (0-K transcripts) to log scaled frequency, where the log_count_bins gives the range for each count bin
#' @param dist a list of master equation simulation results, each element is a vector of length K
#' @param log_count_bins a vector of form seq(min,max,stepsize), or doesn't have equal distance bins
#' @return a matrix where each column is a bin, and each row is one distribution, and the contents are frequencies of probability of being in each bin 
Sim_LogDist <- function(dist,log_count_bins){
    bins=10^log_count_bins
    Log_dist=lapply(dist,function(X){
        inbins=split(X[c(2:length(X))],cut(c(2:length(X)),bins))
        dist=c(X[1],sapply(inbins,sum))
        dist[is.na(dist)]=0
        return(dist)
    })
    Log_dist=do.call(rbind,Log_dist)
    return(Log_dist)
}

#' Getting logged expression distribution
#'
#' Prepares for plotting the Count Heatmap
#' @param dist the expression matrix
#' @param log_count_bins a vector of the cut-offs for the histogram
#' @return a matrix where the rows are the genes and columns are the number of samples within a count category
LogDist <- function(counts,log_count_bins){
    log_dist=apply(log(counts+1,base=10),1,function(x){
        if(sum(is.na(x))!=length(x)){
            dist0=sum(x==0)
            range_c=log_count_bins
            count=table(cut(x[x>0],range_c))
            dist=c(dist0,count)/(sum(count)+dist0)
            return(dist)}else{
                return(NA)
        }
    })
    log_dist=t(log_dist)
    return(log_dist)
}

#' Plotting logged expression distribution
#'
#' takes an expression matrix and makes a 2D histogram on the log scale where each row is a gene and the number of samples in a bin is shown by the intensity of the color
#' @param log_real the logged distribution of count distribution
#' @param mean_counts the average expression for each gene, used for sorting purpose
#' @param given_ord the given order of genes 
#' @param zeropropthres the genes with zeroproportion greater than this number is not plotted (default to 0.8)
#' @param filename the name of the output plot
#' @param saving if the plot should be saved into a file
#' @param data_name a string which is included in the title of the plot to describe the data used
PlotCountHeatmap <- function(log_real, mean_counts,given_ord=NA,zeropropthres=0.8,
                             filename,saving=F, data_name){
    mean_counts=mean_counts[log_real[,1]<zeropropthres]
    log_real=log_real[log_real[,1]<zeropropthres,]
    colnames(log_real)[1]='.0'
    plot_real=melt(log_real)
    plot_real$freq=plot_real$value
    genenames=rownames(log_real)
    if(is.null(genenames)){
        genenames=as.character(c(1:length(log_real[,1])))
        rownames(log_real)=genenames
    }
    if(is.na(given_ord[1])){    
        cluster_ord <- hclust( dist(log_real, method = "euclidean"), method = "ward.D" )$order
        ord=order(cut(log(mean_counts+1),30),order(cluster_ord))
    }else{ord<-given_ord}
    plot_real$Gene <- factor( plot_real$X1, levels = genenames[ord])
    p <- ggplot(plot_real, aes(X2, Gene)) + geom_tile(aes(fill = freq)) +
    scale_fill_gradient(low = "white", high = "black",trans='identity')+
    ggtitle(sprintf('distribution of mRNA counts of %s', data_name)) +
    labs(colour = 'Percentage of Cells',x='log10(Count) bins',y='Genes') +
    scale_y_discrete(breaks=NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if(saving==T){ggsave(filename,dev='jpeg',width = 8, height = 8)}else{return(list(ord,p))}
}
#' Plotting the histograms of kon,koff,s values
#'
#' plot colored histograms of parameters
#' @param params a matrix of 3 columns, the first one is kon, the second is koff and the third is s
#' @param samplename the prefix of the plot, the suffix is '.params_dist.jpeg'
#' @param saving if the plot should be saved to file
#' @return make a plot of three histograms
PlotParamHist<-function(params,samplename,saving=F){
        df <- data.frame(kon = log(base=10,params[,1]),koff=log(base=10,params[,2]),s=log(base=10,params[,3]))
        df <- melt(df)
        p1 <- ggplot(df,aes(x=value)) +
            geom_histogram(data=subset(df,variable == 'kon'),aes(y = ..density..), binwidth=density(df$value)$bw) +
            geom_density(data=subset(df,variable == 'kon'),fill="red", alpha = 0.2) 
        p2 <- ggplot(df,aes(x=value)) +
            geom_histogram(data=subset(df,variable == 'koff'),aes(y = ..density..), binwidth=density(df$value)$bw) +
            geom_density(data=subset(df,variable == 'koff'),fill="green", alpha = 0.2) 
        p3 <- ggplot(df,aes(x=value)) +
            geom_histogram(data=subset(df,variable == 's'),aes(y = ..density..), binwidth=density(df$value)$bw) +
            geom_density(data=subset(df,variable == 's'),fill="blue", alpha = 0.2) 
        if(saving==T){ggsave(paste(samplename,'.params_dist.jpeg',sep=''),plot=arrangeGrob(p1, p2, p3, ncol=1),device='jpeg')}else{p <- arrangeGrob(p1, p2, p3, ncol=1)}
        return(p)
}

#' rescale2range
#'
#' Subfunction for Plotting FNR. Rescale the values in vec such that the lagest is n, and the smallest is 1.
#' @param vec input vector
#' @param n the largest integer to scale the vector to (the smallest is 1)
rescale2range <- function(vec, n){
  a <- (n-1)/(max(vec)-min(vec))
  return(a*vec+(1-a*min(vec)))
}

#' Plotting PCA results (PC1 and PC2)
#' @param PCAres the PCA results
#' @param col_vec a vector to specify the colors for each point
#' @param figuretitle title for the plot
plotPCAbasic <- function(PCAres, col_vec, figuretitle) {
  variance_perc <- 100*(PCAres$sdev)^2/sum((PCAres$sdev)^2)
  plot(PCAres$x[,1], PCAres$x[,2], col=col_vec, pch=20,
       xlab=sprintf("PC1 %4.2f%%", variance_perc[1]), 
       ylab=sprintf("PC2 %4.2f%%", variance_perc[2]),
       main=figuretitle)
}

#' arrange multiple plots from ggplot2 to one figure. 
#' From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' @param cols Number of columns in layout
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



