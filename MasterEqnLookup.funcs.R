##############################################################
#Functions
##############################################################
library('ggplot2')
library('gridExtra')
library("fAsianOptions")
library('FNN')
library('diptest')
# MasterEqn <- function(rateParams,scale=1){
#   s <- rateParams["s"]
#   d <- rateParams["d"]
#   k_on <- rateParams["k_on"]
#   k_off <- rateParams["k_off"]
#   ssMean <- k_on / (k_on + k_off) * s
#   ssVar <- ssMean + ((k_on * k_off)/(k_on + k_off)^2) * s^2 / (k_on + k_off + 1)
#   if(ssVar>1000){ssVar=1000}
#   upperlimit=(ceiling(ssMean + scale*ssVar)+1)
#   # if(upperlimit>3000){upperlimit=3000}
#   x=seq(0, upperlimit, 1)
#   #y <- ((s/d)^x * exp(-s/d) * gamma(k_on/d+x) * gamma(k_on/d+k_off/d)) / 
#   #  (factorial(x) * gamma(k_on/d+k_off/d+x) * gamma(k_on/d)) 
#   # use lgamma and lfactorial to allow values over 170
#   # when s is large, lgamma and lfactorial is still not sufficient (over flow for z?)
#   y <- (x*log(s/d) + (-s/d) + lgamma(k_on/d+x) + lgamma(k_on/d+k_off/d)) -
#     lfactorial(x) - lgamma(k_on/d+k_off/d+x) - lgamma(k_on/d)
#   z <- sapply(x, function(xx)  Re(kummerM(s/d, k_off/d, (k_on/d+k_off/d+xx),lnchf=1)))
#   logp <- y + z
#   p=exp(logp)
#   if (any(is.na(p))){
#     p <- p[-which(is.na(p))]
#   }
#   return(p) # data_1gene is indices
# }
#changes made: 
# 1. instead of rounding ssMean+scale*ssVar, take ceiling (for when there is only 2 categories 0 and 1)
# 2. add one to the range just to be sure
# 3. make kummerM return log scaled values for large s
# 4. cap ssVar at 1000 so that we don't calculate exceedingly small values
# instead of calculating range, just do it at 100 

MasterEqn <- function(rateParams,scale=2){
  k_on <- rateParams[1]
  k_off <- rateParams[2]
  s <- rateParams[3]
  d <- rateParams[4]
 ssMean <- k_on / (k_on + k_off) * s
  ssVar <- ssMean + ((k_on * k_off)/(k_on + k_off)^2) * s^2 / (k_on + k_off + 1)
  if(ssVar>1000){ssVar=1000}
  x=seq(0, (ceiling(ssMean + scale*ssVar)+1), 1)
  #y <- ((s/d)^x * exp(-s/d) * gamma(k_on/d+x) * gamma(k_on/d+k_off/d)) / 
  #  (factorial(x) * gamma(k_on/d+k_off/d+x) * gamma(k_on/d)) 
  # use lgamma and lfactorial to allow values over 170
  # when s is large, lgamma and lfactorial is still not sufficient (over flow for z?)
	y <- (x*log(s/d) + (-s/d) + lgamma(k_on/d+x) + lgamma(k_on/d+k_off/d)) -
    lfactorial(x) - lgamma(k_on/d+k_off/d+x) - lgamma(k_on/d)
	z <- sapply(x, function(xx)  Re(kummerM(s/d, k_off/d, (k_on/d+k_off/d+xx),lnchf=1)))
  logp <- y + z
  p=exp(logp)
  if (any(is.na(p))){
    p <- p[-which(is.na(p))]
  }
  return(p) # data_1gene is indices
}

MasterEqn2 <- function(rateParams,scale=1){
  k_on <- rateParams[1]
  k_off <- rateParams[2]
  s <- rateParams[3]
  d <- rateParams[4]
  p_all=NA
  k=ceiling(s/5)
  for(i in c(0:1000)){
  	if(sum(p_all,na.rm=T)>0.999){break}
	  x=seq(i*k, ((i+1)*k-1), 1)
	  y <- (x*log(s/d) + (-s/d) + lgamma(k_on/d+x) + lgamma(k_on/d+k_off/d)) -
	    lfactorial(x) - lgamma(k_on/d+k_off/d+x) - lgamma(k_on/d)
	  z <- sapply(x, function(xx)  Re(kummerM(s/d, k_off/d, (k_on/d+k_off/d+xx),lnchf=1)))
	  logp <- y + z
	  p=exp(logp)
	  if (any(is.na(p))){
	    p <- p[-which(is.na(p))]
	  }
	  p_all=c(p_all,p)  	
  }
  return(p_all[-1]) # data_1gene is indices
}
DistFeatures <- function(x){
	first=sum(x*c(0:(length(x)-1)))
	second=sum(x*c(0:(length(x)-1))^2)
	third=sum(x*c(0:(length(x)-1))^3)
	sigma2=(second-first^2)^0.5
	skewedness=third/(second)^1.5
	dipp=my.dip.test(x)
	return(c(first,sigma2,skewedness,x[1],dipp))
}
my.dip.test <- function(X){
	delta=diff(X)
	delta[abs(delta)<1e-7]=0
	len=length(X)
	ninfl=0;uptick=0;prop2=0;mode2=0;lastsign=(delta[1]<0)
	if(sum(delta<0)!=(len-1)){
		for(j in c(2:(len-1))){
			# temp[(j-1),]=(c(delta[j-1],lastsign,delta[j],delta[j]<0))
			if( ((delta[j]<0) != lastsign) & delta[j]!=0){
				ninfl=ninfl+1
			}
			if(delta[j]!=0){lastsign = (delta[j]<0)}
		}
		# uptick=min(c(1:(len-1))[delta>0])
		# prop2=sum(X[uptick:len])
		# mode2=c(2:len)[X[2:len]==max(X[2:len])]
	}
	# return(c(ninfl,uptick,prop2,mode2))
	return(ninfl)
}
DistMoments <- function(x,k){
	order=c(1:k)
	i=c(0:(length(x)-1))
	moments=sapply(order,function(Y){sum(x*(i^Y))^(1/k)})
	return(moments)
}

KLdist <- function(a,b){  
	lendiff=length(a)-length(b)
	k=3
	if(lendiff<0){a=c(a,rep(0,abs(lendiff)));k=2}
	if(lendiff>0){b=c(b,rep(0,abs(lendiff)));k=1}
	kl.dist(a,b)[[k]]
}

MasterEqn_heatmap <- function(FUN,data,filename,funcname,transformer,saving=F){
	feature_value=FUN(data)
	plot_data=lapply(unique(data[,3]),function(s_value){
		data.frame(kon=data[data[,3]==s_value,1],koff=data[data[,3]==s_value,2],s=data[data[,3]==s_value,3],value=feature_value[data[,3]==s_value])
	})
	plist = lapply(plot_data,function(X){
		ggplot(X, aes(x=log(kon,base=10),y=log(koff,base=10),fill = value)) + 
		geom_tile() +
		scale_fill_gradientn(colours = rainbow(10),trans = transformer ) + 
		ggtitle(paste('heatmap of ',funcname,', s=',X$s[1],sep='')) +
		labs(colour = funcname)
	})
	if(saving==T){ggsave(filename, marrangeGrob(grobs = plist, nrow=1, ncol=1))}else{
		return(plist)
	}
}

Fano <- function(data){
	fano=data[,5]^2 / data[,4]
	return(fano)
}

FirstMoment <- function(data){
	mean=data[,4]
	return(mean)
}
Sigma <- function(data){
	std=data[,5]
	return(std)
}
Skewedness <- function(data){
	skew=data[,6]
	return(skew)
}
PropZero <- function(data){
	prop=data[,7]
	return(prop)
}
Diptest <- function(data){
	dipp=data[,8]
	return(dipp)
}


