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
#' @examples
#' 
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
#' @param zeropropthres the genes with zeroproportion greater than this number is not plotted (default to 0.8)
#' @param filename the name of the output plot
#' @examples
#' Sim_LogDist()
PlotCountHeatmap <- function(log_real, mean_counts,given_ord=NA,zeropropthres=0.8,filename,saving=F){
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
    ggtitle('distribution of mRNA counts') +
    labs(colour = 'Percentage of Cells',x='Counts',y='Genes') +
    scale_y_discrete(breaks=NULL)
    if(saving==T){ggsave(filename,dev='jpeg',width = 8, height = 8)}else{return(p)}
}
#' Plotting the histograms of kon,koff,s values
#'
#' plot colored histograms of parameters
#' @param params a matrix of 3 columns, the first one is kon, the second is koff and the third is s
#' @param samplename the prefix of the plot, the suffix is '.params_dist.jpeg'
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
