#' Get the best matched parameter 
#'
#' This function matches a real dataset to a database of summary information of simulated datasets, plots a qqplot for user to visualize similarity between their dataset and the simulated dataset, and suggests parameters to use in the simulation
#' @param tech 'ss2','umi' the match database are constructed based on reasonable values for each technology
#' @param counts expression matrix
#' @param plotfilename output name for qqplot
#' @return three set of best matching parameters that was used to simulate the best matching dataset to the experimental dataset

BestMatchParams <- function(tech,counts,plotfilename){
  mean_exprs <- quantile(rowMeans(counts+1,na.rm=T),seq(0,1,0.002))
  fano_exprs <- quantile(apply(counts,1,fano),seq(0,1,0.002),na.rm=T)
  percent0 <- quantile(apply(counts,1,percent_nonzero),seq(0,1,0.002))
  if(tech=='ss2'){
    load('SymSim/grid_summary/exp_figure4ss2_Lgrid.summary.robj')   
  }else if(tech =='umi'){
    load('SymSim/grid_summary/exp_figure4umi_Lgrid.summary.robj')   
  }
  grid_summary <- list(fano_bins,mean_bins,nonzero_bins)

  exp_summary <- list(fano_exprs,mean_exprs,percent0)
  dists <- lapply(c(1:3),function(i){
   dist <- apply(grid_summary[[i]],1,function(X){sum(abs(X-exp_summary[[i]]))})
   return(dist)
  })
  dists <- do.call(cbind,dists)
  dists <- apply(dists,2,function(X){X/mean(X)})
  dists <- rowSums(dists)
  best_match <- c(1:length(dists))[rank(dists)%in% c(1:3)]

  best_params <- lapply(best_match,function(X){sim_params[X,]})
  plotnames <- c('fano','mean','percent_nonzero')
  pdf(file=plotfilename)
  par(mfrow=c(3,3))
  for(i in c(1:3)){
    # for(j in c(1:3)){
      for(k in c(1:3)){
        bin1=grid_summary[[k]][best_match[i],]
        bin2=exp_summary[[k]]
        if(k %in% c(1,2)){bin1 <- log(base=10,bin1);bin2 <- log(base=10,bin2)}
        plot(bin1,bin2,pch=16,xlab='simulated values',ylab='experimental values',main=paste(i,'best','match', plotnames[k]))
        lines(c(-10,10),c(-10,10),col='red')      
      }
    }
  # }
  dev.off()
  best_params <- do.call(rbind,best_params)
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
PlotCountHeatmap <- function(log_real, mean_counts,given_ord=NA,zeropropthres=0.8,filename,saving=F, data_name){
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
#' Calculating the 
#'
#' @param filename the name of the simualted robj data file

cv <- function(x) {return(sd(x)/mean(x))}
fano <- function(x) {return((sd(x))^2/mean(x))}

percent_nonzero <- function(x) {return(sum(x>0)/length(x))}



#' Plotting simulated FNR 
#'
#' @param current_counts expression matrix
#' @param true_counts true expression matrix
#' @param titlestr title string

plotFNRsim <- function(current_counts, true_counts, x_use="true", titlestr){
  
  #cellsize <- colSums(current_counts)
  zeros <- colSums(current_counts==0 & true_counts!=0)
  #summary(zeros)
  
  ref.glms <- lapply(c(1:length(true_counts[1,])),function(i){
    real <- true_counts[,i]
    bias <- current_counts[,i]/sum(current_counts[,i])*10^6
    fn <- (bias==0 & real!=0)
    if (x_use == "observed")
    {data <- data.frame(exprs=log(bias+1,10), fn=fn)} else
      if (x_use == "true")
      {data <- data.frame(exprs=log(real+1,10), fn=fn)}
    data <- data[real!=0,]
    model <- glm(fn~exprs,family=binomial(link='logit'),data=data)
    return(model$coefficients)
  })
  
  ncols <- 2000
  colfunc <- colorRampPalette(c("blue", "green", "red", "purple"))
  allcol <- colfunc(ncols)
  
  col_vec <- allcol[round(rescale2range(colSums(current_counts), 2000))]
  
  plot(NULL, main = sprintf("FNR Curves %s", titlestr), ylim = c(0,1),xlim = c(0,6), 
       ylab = "Failure Probability", xlab = "Mean log10 Expression")
  x = (0:60)/10
  for(si in 1:length(ref.glms)){
    y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
    lines(x, y, type = 'l', lwd = 2, col=col_vec[si])
  }
}

#' Plotting real data FNR 
#'
#' @param expr_matrix expression matrix
#' @param hk house keeping genes
#' @param data_name title string
#' @param col_vec color vector

plotFNRreal <- function(expr_matrix, hk, data_name, col_vec){
  # Mean log10(x+1) expression
  mu_obs = rowMeans(log10(expr_matrix[hk,]+1))
  drop_outs = expr_matrix[hk,] == 0
  
  # Logistic Regression Model of Failure
  ref.glms = list()
  for (si in 1:dim(drop_outs)[2]){
    #fit = glm(cbind(drop_outs[,si], 1 - drop_outs[,si]) ~ mu_obs,family=binomial(logit))
    fit = glm(drop_outs[,si] ~ mu_obs,family=binomial(logit))
    ref.glms[[si]] = fit$coefficients
  }
  
  #The list ref.glm contains the intercept and slope of each fit. 
  # We can now visualize the fit curves and the corresponding Area Under the Curves (AUCs):
  plot(NULL, main = sprintf("FNR Curves %s", data_name), ylim = c(0,1),xlim = c(0,6), 
       ylab = "Failure Probability", xlab = "Mean log10 Expression")
  x = (0:60)/10
  for(si in 1:ncol(expr_matrix)){
    y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
    lines(x, 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1), type = 'l', lwd = 2, col=col_vec[si])
  }
}


plotPCAbasic <- function(PCAres, col_vec, figuretitle) {
  variance_perc <- 100*(PCAres$sdev)^2/sum((PCAres$sdev)^2)
  plot(PCAres$x[,1], PCAres$x[,2], col=col_vec, pch=20,
       xlab=sprintf("PC1 %4.2f%%", variance_perc[1]), 
       ylab=sprintf("PC2 %4.2f%%", variance_perc[2]),
       main=figuretitle)
}

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

ks.dist<- function (spec1, spec2, f = NULL, mel = FALSE, plot = FALSE, 
    type = "l", lty = c(1, 2), col = c(2, 4), flab = NULL, alab = "Cumulated amplitude", 
    flim = NULL, alim = NULL, title = TRUE, legend = TRUE) 
{
    leg <- c(as.character(deparse(substitute(spec1))), as.character(deparse(substitute(spec2))))
    if (is.matrix(spec1) && ncol(spec1) == 2) 
        s1 <- spec1[, 2]
    else s1 <- spec1
    if (is.matrix(spec2) && ncol(spec2) == 2) 
        s2 <- spec2[, 2]
    else s2 <- spec2
    n1 <- length(s1)
    n2 <- length(s2)
    if (n1 != n2) 
        stop("spec1 and spec2 must have the same length")
    if (any(is.na(s1)) | any(is.na(s2))) {
        D <- F <- NA
        warning("The data set contains 'NA' values. The returned values have been set to NA.", 
            call. = FALSE)
        res <- list(D = D, F = F)
        return(res)
    }
    else {
        if (any(s1 < 0) | any(s2 < 0)) 
            stop("spectra (spec 1 and/or spec 2) do not have to be in dB")
        if (sum(s1) == 0) {
            warning(paste("Caution!, spec1 is a null spectrum"), 
                call. = FALSE)
        }
        if (sum(s2) == 0) {
            warning(paste("Caution!, spec2 is a null spectrum"), 
                call. = FALSE)
        }
        s1 <- s1/sum(s1)
        s2 <- s2/sum(s2)
        cum.s1 <- cumsum(s1)
        cum.s2 <- cumsum(s2)
        diff <- abs(cum.s1 - cum.s2)
        D <- max(diff)
        if (!is.null(f) & mel) {
            f <- 2 * mel(f/2)
        }
        if (is.null(f)) {
            if (is.vector(spec1) & is.vector(spec2)) {
                if (plot) 
                  warning("'f' is missing, F cannot be computed - No plot produced", 
                    call. = FALSE)
                else warning("'f' is missing, F cannot be computed", 
                  call. = FALSE)
                res <- list(D = D, F = NA)
                return(res)
            }
            else {
                if (is.matrix(spec1)) 
                  f <- spec1[nrow(spec1), 1] * 2000 * nrow(spec1)/(nrow(spec1) - 
                    1)
                else if (is.matrix(spec2)) 
                  f <- spec2[nrow(spec2), 1] * 2000 * nrow(spec2)/(nrow(spec2) - 
                    1)
            }
        }
        x <- seq(0, f/2 * (n1 - 1)/n1, length.out = n1)/1000
        pos <- which.max(diff)
        F <- x[pos]
        res <- list(D = D, F = F)
        if (plot) {
            if (mel) {
                w <- " kmel"
            }
            else {
                w <- " kHz"
            }
            if (is.null(alim)) 
                alim <- c(0, 1)
            if (is.null(flim)) 
                flim <- c(0, f/2000)
            if (is.null(flab)) {
                if (mel) 
                  flab <- "Frequency (kmel)"
                else flab <- "Frequency (kHz)"
            }
            x <- seq(0, (f/2) * (n1 - 1)/n1, length.out = n1)/1000
            plot(x = x, y = cum.s1, col = col[1], lty = lty[1], 
                type = type, xaxs = "i", xlab = flab, xlim = flim, 
                yaxs = "i", ylim = alim, ylab = alab, ...)
            lines(x = x, y = cum.s2, col = col[2], lty = lty[2], 
                type = type)
            segments(x0 = F, y0 = cum.s1[pos], x1 = F, y1 = cum.s2[pos], 
                lwd = 2)
            if (title) 
                title(main = paste("D = ", round(D, 3), "\\n F = ", 
                  round(F, 3), w))
            if (legend) 
                legend("topleft", col = col, lty = lty, legend = leg, 
                  bty = "n")
            invisible(res)
        }
        else return(res)
    }
}
