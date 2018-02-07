library('devtools')
load_all("SymSim")
dir <- '/data/yosef2/users/xiuwei/simulation/Rvariables/exp_3pop0205/'

files <- list.files(dir)
k <- sapply(strsplit(files,'_'),function(X){
	X <- strsplit(X[2],'.',fixed=T)[[1]][1]
	as.numeric(X)
})
files <- files[order(k)]

pop1=1;pop2=3;check_params=T

lapply(files,function(filename){
	load(paste(dir,filename,sep=''))
	true_counts <- true_counts_res[[1]]
	GE <- true_counts_res[[2]]
	meta <- true_counts_res[[3]]
	params <- true_counts_res[[4]]

	temp <- colnames(meta)
	DE_evfs <- lapply(c('kon','koff','s'),function(i){
		k <- c(1:length(temp))[grep(i,temp)]
		nonDE <- k[grep('non_DE',temp[k])]
		# this is temporary for datasets simulated with mistakenly marked EVFs
		DE <- k[k%in%nonDE]
		DE <- DE+1-min(k)
		return(DE)
	})

	nonDEgenes <- lapply(c(1:3),function(parami){
		id <- rowSums(GE[[parami]][,DE_evfs[[parami]]]!=0) ==0
		if(check_params ==T){
			de_genes <- params[[parami]][id,]				
		}else{de_genes <- true_counts[id,]}
		return(de_genes)
	})

	DEgenes <- lapply(c(1:3),function(parami){
		id <- rowSums(GE[[parami]][,DE_evfs[[parami]]]!=0) >= floor(length(DE_evfs[[parami]])/2)
		if(check_params ==T){
			de_genes <- params[[parami]][id,]				
		}else{de_genes <- true_counts[id,]}
		return(de_genes)
	})

	DE_GE <- lapply(c(0:min(sapply(DE_evfs,length))),function(k){
		lapply(c(1:3),function(parami){
			id <- rowSums(GE[[parami]][,DE_evfs[[parami]]]!=0) ==k
			if(check_params ==T){
				de_genes <- params[[parami]][id,]				
			}else{de_genes <- true_counts[id,]}
			return(de_genes)
		})
	})

	mean_FC <- lapply(c(1:3),function(parami){
		FC <- sapply(DE_GE,function(X){
			if(is.null(X[[parami]])){return(NA)
			}else if(length(dim(X[[parami]]))==0){return(NA)
			}else{
				FC <- log(rowMeans(X[[parami]][,meta[,2]==pop1])/rowMeans(X[[parami]][,meta[,2]==pop2]))
				return(FC)				
			}
		})
		mean_FC <- sapply(FC,function(X){mean(abs(X[is.finite(X)]))})
		return(mean_FC)
	})

	pdf(file=paste('DE_plots/',filename,'.DE.params.pdf',sep=''))
	par(mfrow=c(2,3))
	for( parami in c(1:3)){
		nonDE_FC <- log(rowMeans(nonDEgenes[[parami]][,meta[,2]==pop1])/rowMeans(nonDEgenes[[parami]][,meta[,2]==pop2]))
		# nonDE_mean <- rowMeans(nonDEgenes[[parami]])
		nonDE_p.value <- sapply(c(1:length(nonDEgenes[[parami]][,1])),function(i){
		wilcox.test(nonDEgenes[[parami]][i,meta[,2]==pop1],nonDEgenes[[parami]][i,meta[,2]==pop2])$p.value
		})
		DE_FC <- log(rowMeans(DEgenes[[parami]][,meta[,2]==pop1])/rowMeans(DEgenes[[parami]][,meta[,2]==pop2]))
		# DE_mean <- rowMeans(DEgenes[[parami]])
		DE_p.value <- sapply(c(1:length(DEgenes[[parami]][,1])),function(i){
		wilcox.test(DEgenes[[parami]][i,meta[,2]==pop1],DEgenes[[parami]][i,meta[,2]==pop2])$p.value
		})	
		plot(x=log(base=10,DE_p.value),y=(DE_FC),pch=16,col='red',main=c('kon','koff','s')[parami])
		points(x=log(base=10,nonDE_p.value),y=(nonDE_FC),pch=16,col='blue')
	}
	for(i in c(1:3)){
		plot(mean_FC[[i]],x=c(0:length(DE_evfs[[i]])),pch=16,type='o',main=c('kon','koff','s')[i])
	}
	dev.off()
})

check_params=F

lapply(files,function(filename){
	load(paste(dir,filename,sep=''))
	true_counts <- true_counts_res[[1]]
	GE <- true_counts_res[[2]]
	meta <- true_counts_res[[3]]
	params <- true_counts_res[[4]]

	temp <- colnames(meta)
	DE_evfs <- lapply(c('kon','koff','s'),function(i){
		k <- c(1:length(temp))[grep(i,temp)]
		nonDE <- k[grep('non_DE',temp[k])]
		DE <- k[!k%in%nonDE]
		DE <- DE+1-min(k)
		return(DE)
	})

	nonDEgenes <- lapply(c(1:3),function(parami){
		id <- rowSums(GE[[parami]][,DE_evfs[[parami]]]!=0) ==0
		if(check_params ==T){
			de_genes <- params[[parami]][id,]				
		}else{de_genes <- true_counts[id,]}
		return(de_genes)
	})

	DEgenes <- lapply(c(1:3),function(parami){
		id <- rowSums(GE[[parami]][,DE_evfs[[parami]]]!=0) >= floor(length(DE_evfs[[parami]])/2)
		if(check_params ==T){
			de_genes <- params[[parami]][id,]				
		}else{de_genes <- true_counts[id,]}
		return(de_genes)
	})

	DE_GE <- lapply(c(0:min(sapply(DE_evfs,length))),function(k){
		lapply(c(1:3),function(parami){
			id <- rowSums(GE[[parami]][,DE_evfs[[parami]]]!=0) ==k
			if(check_params ==T){
				de_genes <- params[[parami]][id,]				
			}else{de_genes <- true_counts[id,]}
			return(de_genes)
		})
	})

	mean_FC <- lapply(c(1:3),function(parami){
		FC <- sapply(DE_GE,function(X){
			if(is.null(X[[parami]])){return(NA)
			}else if(length(dim(X[[parami]]))==0){return(NA)
			}else{
				FC <- log(rowMeans(X[[parami]][,meta[,2]==pop1])/rowMeans(X[[parami]][,meta[,2]==pop2]))
				return(FC)				
			}
		})
		mean_FC <- sapply(FC,function(X){mean(abs(X[is.finite(X)]))})
		return(mean_FC)
	})

	pdf(file=paste('DE_plots/',filename,'.DE.counts.pdf',sep=''))
	par(mfrow=c(2,3))
	for( parami in c(1:3)){
		nonDE_FC <- log(rowMeans(nonDEgenes[[parami]][,meta[,2]==pop1])/rowMeans(nonDEgenes[[parami]][,meta[,2]==pop2]))
		# nonDE_mean <- rowMeans(nonDEgenes[[parami]])
		nonDE_p.value <- sapply(c(1:length(nonDEgenes[[parami]][,1])),function(i){
		wilcox.test(nonDEgenes[[parami]][i,meta[,2]==pop1],nonDEgenes[[parami]][i,meta[,2]==pop2])$p.value
		})
		DE_FC <- log(rowMeans(DEgenes[[parami]][,meta[,2]==pop1])/rowMeans(DEgenes[[parami]][,meta[,2]==pop2]))
		# DE_mean <- rowMeans(DEgenes[[parami]])
		DE_p.value <- sapply(c(1:length(DEgenes[[parami]][,1])),function(i){
		wilcox.test(DEgenes[[parami]][i,meta[,2]==pop1],DEgenes[[parami]][i,meta[,2]==pop2])$p.value
		})	
		plot(x=log(base=10,DE_p.value),y=(DE_FC),pch=16,col='red',main=c('kon','koff','s')[parami])
		points(x=log(base=10,nonDE_p.value),y=(nonDE_FC),pch=16,col='blue')
	}
	for(i in c(1:3)){
		plot(mean_FC[[i]],x=c(0:length(DE_evfs[[i]])),pch=16,type='o',main=c('kon','koff','s')[i])
	}
	dev.off()
})
