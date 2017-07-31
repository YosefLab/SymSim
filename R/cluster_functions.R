#' Perform tSNE 
#' 
#' This function takes the ground truth (meta), the expresssion matrix (data), and performs tSNE, saves the result to a jpeg file with user specified name
#' @param meta: simulation parameters
#' @param data: expression matrix
#' @param plotname: the name of the jpeg file
#' @param label: the column name of the meta data that the points needs to be colored by

PlotTsne <- function(meta,data,plotname,label,discrete=T,saving=F){
	library('Rtsne')
    uniqcols<-c(1:length(data[1,]))[!duplicated(t(data))]
    data <- data[,uniqcols];meta <- meta[uniqcols,] 
    uniqrows<-c(1:length(data[,1]))[!duplicated(data)]
    data <- data[uniqrows,]
	data_tsne=Rtsne(t(data))
	if(discrete==T){
		plot_tsne=data.frame(alpha=meta[,'alpha'],beta=meta[,'beta'],sigma=meta[,'sigma'],
		pop=meta[,'pop'],label=factor(meta[,label]),x=data_tsne$Y[,1],y=data_tsne$Y[,2])		
	}else{
		plot_tsne=data.frame(alpha=meta[,'alpha'],beta=meta[,'beta'],sigma=meta[,'sigma'],
		pop=meta[,'pop'],label=meta[,label],x=data_tsne$Y[,1],y=data_tsne$Y[,2])
	}
	p <- ggplot(plot_tsne, aes(x, y))
	p <- p + geom_point()
	p <- p + geom_point(aes(colour = plot_tsne[['label']]))+labs(color=label)
	if(saving==T){ggsave(p,filename=plotname,device='jpeg',width=5,height=4)}else{print(p)}
	return(plot_tsne)
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
		plot_pca=data.frame(alpha=meta[,'alpha'],beta=meta[,'beta'],sigma=meta[,'sigma'],
		pop=meta[,'pop'],label=factor(meta[,label]),x=data_pc$x[,1],y=data_pc$x[,2])		
	}else{
		plot_pca=data.frame(alpha=meta[,'alpha'],beta=meta[,'beta'],sigma=meta[,'sigma'],
		pop=meta[,'pop'],label=meta[,label],x=data_pc$x[,1],y=data_pc$x[,2])
	}
	p <- ggplot(plot_pca, aes(x, y))
	p <- p + geom_point()
	p <- p + geom_point(aes(colour = plot_pca[['label']]))+labs(color=label)
	if(saving==T){ggsave(p,filename=plotname,device='jpeg',width=5,height=4)}else{print(p)}
	return(plot_pca)
}
