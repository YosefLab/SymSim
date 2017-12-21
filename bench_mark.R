#######################
# Clustering 
#######################

required_packages <- c('devtools','FNN','reshape','ggplot2','gridExtra','ape','MASS','Rtsne','RColorBrewer',
                       'fOptions','fBasics','timeSeries','timeDate','Biobase','repr' )
install_packages <- rownames(installed.packages())
need_install <- setdiff(required_packages,install_packages)
if(length(need_install)>0){install.packages(need_install)}

load_result <- lapply(required_packages, require, character.only = TRUE)

load('allsim.20170706.robj')
load('param_realdata.robj')
load_all('SCsimR')



param_names <- c('bimod','alpha','evf_sd','batch','robj')


files <- list.files('20170804_sim/')
files <- files[grep('robj',files)]

param_values <- lapply(c(1:(length(param_names)-1)),function(name){
	param_value <- sapply(files,function(X){
		X=strsplit(X,paste('.',param_names[name],'=',sep=''),fixed=T)[[1]][2]
		X=strsplit(X,paste('.',param_names[name+1],sep=''),fixed=T)[[1]][1]
		return(as.numeric(X))
	})	
	return(param_value)
})
param_values = do.call(cbind,param_values)
rownames(param_values)=NULL
colnames(param_values)=param_names[c(1:4)]

files = files[param_values[,4]==0.5 & param_values[,3]==0.2]

#######################
# calculate performance for 1 batch files 
#######################
files <- files[grep('1batch',files)]
sink('benchmark.0804.txt')
library(fossil)
print('start')
clustering_perf <- lapply(files,function(X){
	load(paste('20170804_sim/',X,sep=''))
	# ptm <- proc.time()
	tsne <- PlotTsne(meta=result[[2]],data=result[[1]][[3]][c(1:1000),],label='pop',plotname='')
	estimated <- as.numeric(kmeans(tsne[[1]][,c(6,7)],5)$cluster)
	true <- as.numeric(tsne[[1]]$pop)
	all_batch <- c(rand.index(true,estimated),adj.rand.index(true,estimated))	
	randi <- rbind(do.call(rbind,per_batch),all_batch)
	# proc.time() - ptm
	print(randi)
	return(randi)
})
sink()

#######################
# calculate performance for 3 batches files 
#######################
sink('benchmark.txt')
library(fossil)
print('start')
clustering_perf <- lapply(files,function(X){
	load(paste('20170819_sim/',X,sep=''))
	# ptm <- proc.time()
	per_batch <- lapply(c(1:3),function(batch){
		meta <- result[[2]]
		tsne <- PlotTsne(meta=meta[meta$batch==batch,],data=result[[1]][[3]][c(1:1000),meta$batch==batch],label='pop',plotname='')
		estimated <- as.numeric(kmeans(tsne[[1]][,c(6,7)],5)$cluster)
		true <- as.numeric(tsne[[1]]$pop)
		c(rand.index(true,estimated),adj.rand.index(true,estimated))
		# print(proc.time() - ptm)
	})
	tsne <- PlotTsne(meta=result[[2]],data=result[[1]][[3]][c(1:1000),],label='pop',plotname='')
	estimated <- as.numeric(kmeans(tsne[[1]][,c(6,7)],5)$cluster)
	true <- as.numeric(tsne[[1]]$pop)
	all_batch <- c(rand.index(true,estimated),adj.rand.index(true,estimated))	
	randi <- rbind(do.call(rbind,per_batch),all_batch)
	# proc.time() - ptm
	print(randi)
	return(randi)
})
sink()
plot(x=0,y=0,xlim=c(0,1),ylim=c(0.2,0.75),pch='.',xlab='EVF variance',ylab='adjusted random index')
col <- rainbow(length(unique(param_values[,2])))
coli <- 1
for (bimod in unique(param_values[,1])){
	for (alpha in unique(param_values[,2])){
		coli <- coli+1
		for (batch in unique(param_values[,4])){
			exp <- (param_values[,1]==bimod & 
				param_values[,2]==alpha & 
				param_values[,4]==batch)
			x <-param_values[exp,3]
			y <- sapply(clustering_perf[exp],function(X){mean(X[c(1:3),2])})
			points(x=x,y=y,type='l',col=col[coli],lwd=2)
		}
	}
}
#######################
# plotting multiple line trend plot
#######################
plot(x=0,y=0,xlim=c(0,1),ylim=c(0.1,0.75),pch='.',xlab='Batch Effect',ylab='adjusted random index')
col <- rainbow(length(unique(param_values[,2])))
coli <- 1
for (bimod in unique(param_values[,1])){
	for (alpha in unique(param_values[,2])){
			coli <- coli+1
		for (sigma in unique(param_values[,3])){
			exp <- (param_values[,1]==bimod & 
				param_values[,2]==alpha & 
				param_values[,3]==sigma)
			x <-param_values[exp,4]
			y <- sapply(clustering_perf[exp],function(X){X[4,2]})
			points(x=x,y=y,type='l',col=col[coli],lwd=2)
		}
	}
}

#######################
# Trajectory 
#######################


