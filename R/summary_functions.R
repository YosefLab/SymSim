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

PlotCountHeatmap <- function(log_real, mean_counts,zeropropthres=0.8,filename){
    mean_counts=mean_counts[log_real[,1]<zeropropthres]
    log_real=log_real[log_real[,1]<zeropropthres,]
    cluster_ord <- hclust( dist(log_real, method = "euclidean"), method = "ward.D" )$order
    genenames=rownames(log_real)
    if(is.null(genenames)){
        genenames=as.character(c(1:length(log_real[,1])))
        rownames(log_real)=genenames
    }
    colnames(log_real)[1]='0'
    plot_real=melt(log_real)
    plot_real$freq=plot_real$value
    ord=order(cut(log(mean_counts+1),30),order(cluster_ord))
    plot_real$Gene <- factor( plot_real$Var1, levels = genenames[ord])
    p <- ggplot(plot_real, aes(Var2, Gene)) + geom_tile(aes(fill = freq)) +
        scale_fill_gradient(low = "white", high = "black",trans='identity')+ 
        ggtitle('distribution of mRNA counts') +
        labs(colour = 'Percentage of Cells',x='Counts',y='Genes') + 
        scale_y_discrete(breaks=NULL)
    ggsave(filename,dev='jpeg',width = 8, height = 8)
}