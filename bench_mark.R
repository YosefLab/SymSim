#######################
# Clustering 
#######################

<<<<<<< 176a52e550dca83672076caa19801f69fc1b3d4b
<<<<<<< HEAD
=======
>>>>>>> update vignette files
required_packages <- c('devtools','FNN','reshape','ggplot2','gridExtra','ape','MASS','Rtsne','RColorBrewer',
                       'fOptions','fBasics','timeSeries','timeDate','Biobase','repr' )
install_packages <- rownames(installed.packages())
need_install <- setdiff(required_packages,install_packages)
if(length(need_install)>0){install.packages(need_install)}

load_result <- lapply(required_packages, require, character.only = TRUE)

load('allsim.20170706.robj')
load('param_realdata.robj')
<<<<<<< 176a52e550dca83672076caa19801f69fc1b3d4b
load_all('SCsimR')



param_names <- c('bimod','alpha','evf_sd','batch','robj')
=======
library('devtools')
# load('allsim.20170706.robj')
# load('match_params.robj')

load_all('SCsimR')
# param_names <- c('bimod','alpha','evf_sd','batch','3batches')
# param_names <- c('bimod','alpha','evf_sd','batch','robj')
param_names <- c('beta','alpha','evf_sd','batch','discrete')
>>>>>>> 99929bdf19520594c788b04b39fb5c4476efbe77
=======
load_all('SCsimR')



param_names <- c('bimod','alpha','evf_sd','batch','robj')
>>>>>>> update vignette files


files <- list.files('20170905_sim/')
files <- files[grep('robj',files)]
files <- files[grep('discrete',files)]
files <- files[grep('20170905.beta',files)]

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


#######################
# calculate performance for 1 batch files 
#######################
param_values <- param_values[grep('1batch',files),]
files <- files[grep('1batch',files)]
files = files[param_values[,4]==0.5 & param_values[,3]==0.2]
param_values <- param_values[param_values[,4]==0.5 & param_values[,3]==0.2,]

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
	print(all_batch)
	return(all_batch)
})
sink()

jpeg('change_alpha.1batch.jpg')
plot(x=0,y=0,xlim=c(0,0.5),ylim=c(0.5,0.8),pch='.',xlab='alpha',ylab='adjusted random index',cex=5)
col <- rainbow(length(unique(param_values[,1])))
coli <- 1
for (bimod in unique(param_values[,1])){
	coli <- coli+1
	for (batch in unique(param_values[,4])){
		for (sigma in unique(param_values[,3])){
			exp <- (param_values[,1]==bimod & 
				param_values[,4]==batch & 
				param_values[,3]==sigma)
			x <-param_values[exp,2]
			y <- sapply(clustering_perf[exp],function(X){X[2]})
			points(x=x,y=y,type='l',col=col[coli],lwd=2)
		}
	}
}
dev.off()

jpeg('change_bimod.1batch.jpg')
plot(x=0,y=0,xlim=c(0,0.5),ylim=c(0.5,.8),pch='.',xlab='bimod',ylab='adjusted random index',cex=5)
col <- rainbow(length(unique(param_values[,1])))
coli <- 1
for (alpha in unique(param_values[,2])){
	coli <- coli+1
	for (batch in unique(param_values[,4])){
		for (sigma in unique(param_values[,3])){
			exp <- (param_values[,2]==alpha & 
				param_values[,4]==batch & 
				param_values[,3]==sigma)
			x <-param_values[exp,1]
			y <- sapply(clustering_perf[exp],function(X){X[2]})
			y <- y[order(x)]
			x <- x[order(x)]
			points(x=x,y=y,type='l',col=col[coli],lwd=2)
		}
	}
}
dev.off()


#######################
# calculate performance for 3 batches files 
#######################
sink('benchmark.20170830.txt')
library(fossil)
print('start')
clustering_perf <- lapply(files,function(X){
	load(paste('20170830_sim/',X,sep=''))
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
plot(x=0,y=0,xlim=c(0.5,1),ylim=c(0.5,0.9),pch='.',xlab='beta',ylab='adjusted random index')
col <- rainbow(length(unique(param_values[,2])))
coli <- 1

# for (bimod in unique(param_values[,1])){
	for (alpha in unique(param_values[,2])){
		coli <- coli+1
		# for (batch in unique(param_values[,4])){
			# exp <- (param_values[,1]==bimod & 
			# 	param_values[,2]==alpha & 
			# 	param_values[,4]==batch)
			exp <- param_values[,2]==alpha
			x <-param_values[exp,1]
			y <- sapply(clustering_perf[exp],function(X){mean(X[c(1:3),2])})
			points(x=x,y=y,type='l',col=col[coli],lwd=2)
		# }
	# }
}

plot(x=0,y=0,xlim=c(0.5,1),ylim=c(0.5,0.9),pch='.',xlab='beta',ylab='adjusted random index')
col <- rainbow(length(unique(param_values[,1])))

coli <- 1
plot(x=0,y=0,xlim=c(0,0.2),ylim=c(0.5,0.9),pch='.',xlab='alpha',ylab='adjusted random index')
for (bimod in unique(param_values[,1])){
	# for (alpha in unique(param_values[,2])){
		coli <- coli+1
		# for (batch in unique(param_values[,4])){
			# exp <- (param_values[,1]==bimod & 
			# 	param_values[,2]==alpha & 
			# 	param_values[,4]==batch)
			exp <- param_values[,1]==bimod
			x <-param_values[exp,2]
			y <- sapply(clustering_perf[exp],function(X){mean(X[c(1:3),2])})
			y <- y[order(x)]
			x <- x[order(x)]
			points(x=x,y=y,type='l',col=col[coli],lwd=2)
		# }
	# }
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


