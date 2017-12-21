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
    labs(colour = 'Percentage of Cells',x='log10(Count) bins',y='Genes') +
    scale_y_discrete(breaks=NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) ##
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

#' rescale2range
#'
#' Subfunction for Plotting FNR
rescale2range <- function(vec, n){ # rescale the values in vec such that the lagest is n, and the smallest is 1.
  a <- (n-1)/(max(vec)-min(vec))
  return(a*vec+(1-a*min(vec)))
}

#' Plotting FNR 
#'
#' Make 2 plots: one is the FNR curves, and the other one is the barplots of AUC
#' @param expr_matrix expression matrix
#' @param data_name what to name the plot
#' @param ncols number of colors in the color ramp

plotFNR <- function(expr_matrix, data_name, ncols){
    rescale2range <- function(vec, n){ # rescale the values in vec such that the lagest is n, and the smallest is 1.
      a <- (n-1)/(max(vec)-min(vec))
      return(a*vec+(1-a*min(vec)))
    }
    # filter out cells where less than 10% of all genes are expressed
    notzero <- apply(expr_matrix,2,function(X){sum(X>0)>length(X)*0.1})
    expr_matrix<-expr_matrix[,notzero]

    # home keeping gene is defined as genes that are in the top half highest expressed genes
    # taht are also in the top half lowest variance genes
    highexprs <- c(1:length(expr_matrix[,1]))[rowMeans(expr_matrix)>quantile(rowMeans(expr_matrix),0.5)]
    sd2 <- apply(expr_matrix,1,var)
    lowvar <- c(1:length(expr_matrix[,1]))[sd2/rowMeans(expr_matrix)<quantile(sd2/rowMeans(expr_matrix),0.5,na.rm=T)]
    hk <- intersect(highexprs,lowvar)

    # Mean log10(x+1) expression
    colfunc <- colorRampPalette(c("blue", "green", "red", "purple"))
    allcol <- colfunc(ncols)

    col_vec <- allcol[round(rescale2range(colSums(expr_matrix), 2000))]
    mu_obs = rowMeans(log10(expr_matrix[hk,]+1))
    drop_outs = expr_matrix[hk,] == 0

    # Logistic Regression Model of Failure
    ref.glms = list()
    for (si in 1:dim(drop_outs)[2]){
    fit = glm(cbind(drop_outs[,si],1 - drop_outs[,si]) ~ mu_obs,family=binomial(logit))
    ref.glms[[si]] = fit$coefficients
    }
    #The list ref.glm contains the intercept and slope of each fit. 
    # We can now visualize the fit curves and the corresponding Area Under the Curves (AUCs):

    plot(NULL, main = sprintf("FNR Curves %s", data_name), ylim = c(0,1),xlim = c(0,2), 
       ylab = "Failure Probability", xlab = "Mean log10 Expression")
    x = (0:60)/10
    AUC = NULL
    for(si in 1:ncol(expr_matrix)){
    y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
    AUC[si] = sum(y)/10
    lines(x, 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1), type = 'l', lwd = 2, col=col_vec[si])
}
# Barplot of FNR AUC

  o = order(AUC)
  barplot(AUC[o], col=col_vec[o], border=col_vec[o], main="FNR AUC")
  
}


#' Plotting GC bias and Length 
#'
#' @param filename the name of the simualted robj data file
PlotGCLENbias <- function(filename,outputname){
    load(filename)
    GCbias <- result[[4]][[1]][[1]][[1]]
    LNbias <- 1/result[[4]][[1]][[1]][[2]]
    plots <- lapply(list(GCbias,LNbias),function(bias){
        bybins <- lapply(c(1:10),function(i){
            X<-result[[1]][[3]][as.numeric(factor(bias))==i,]
            average<-rowMeans(X)
            average <- average[average>0]
            mu <- mean(log(average,base=2),na.rm=T)
            sd <- sd(log(average,base=2),na.rm=T)
            c(mu,sd,i)
        })
        bybins <- do.call(rbind,bybins)
        colnames(bybins)<-c('mu','sd','bin')
        temp <- as.data.frame(bybins)
        pd <- position_dodge(0.1)
        p1 <- ggplot(temp, aes(x=bin, y=mu)) + 
            geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd), width=.1, position=pd) +
            geom_line(position=pd) +
            geom_point(position=pd)
        p2 <- ggplot(temp,aes(x=bin,y=sd/mu))+
            geom_line(position=pd) +
            geom_point(position=pd)
        return(list(p1,p2))
    })
    plots<-do.call(c,plots)
    names(plots) <- c('GC_mean','GC_fano','Len_mean','Len_fano')
    ggsave(outputname, arrangeGrob(grobs = plots))
    return(plots)
}


