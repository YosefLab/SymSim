#' Calculating Master Equation
#'
#' This function calculates the distribution of transcritp counts in single cells given the dynamic parameters for gene expression
#' @param rateParams a vector of length 4 including kon, koff, s, d
#' @return a vector of X where the probability of having less or equal to X transcripts is greater than 0.999
#' @examples 
#' MasterEqn2(c(1,1,10,0.1))
MasterEqn2 <- function(rateParams){
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

#' Getting True Counts from EVF and Gene effects
#'
#' This function first calls the Get_params function to calculate the kinetic parameters from evf and gene effects. 
#' And then tweaks the bimodality variable
#' At last finds the best matched simulation and sample from the master equation and outputs the true expression of a single cell 
#' @param allparams the sets of parameter used to simulate master equations
#' @param sim_master the master equation results
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param bimdod the parameter for increasing the bimodality of dataset. Takes value between 0 and 1, where 0 means keep the kinetic parameters as is, and 1 means move the kinetic parameters on the kon-koff plane to a space where the distribution of transcript counts is bimodal 
#' @return the simulated (true) expression matrix
#' @examples 
#' EVF2TrueCounts()
EVF2TrueCounts <- function(allparams,matched_params,sim_master,evf,gene_effects,bimod){
	nevf <- length(evf)
	ngenes <- length(gene_effects[[1]][,1])
	params <- Get_params(gene_effects,evf,match_params,bimod)
	best_matches <- FNN::knnx.index(data=allparams,query=params,k=1,algorithm=c('kd_tree'))
	sim <- sim_master[best_matches]
	sim_exprs <- SampleExprs(sim,1)
	return(sim_exprs)
}


#' Getting GeneEffects matrices 
#'
#' This function randomly generates the effect size of each evf on the dynamic expression parameters
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param randomseed (should produce same result if ngenes, nevf and randseed are all the same)
#' @param prob the probability that the effect size is not 0
#' @param sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from 
#' @return a list of 3 matrices, each of dimension ngenes * nevf
#' @examples 
#' GeneEffects()
GeneEffects <- function(ngenes,nevf,randseed,prob,sd){
	set.seed(randseed)
	lapply(c('kon','koff','s'),function(param){
		effect <- lapply(c(1:ngenes),function(i){
			nonzero <- sample(size=nevf,x=c(0,1),prob=c((1-prob),prob),replace=T)
			nonzero[nonzero!=0]=rnorm(sum(nonzero),mean=0,sd=sd)
			return(nonzero)
		})
		return(do.call(rbind,effect))
	})
}

#' Getting GeneEffects matrices 
#'
#' This function randomly generates the effect size of each evf on the dynamic expression parameters
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param randomseed (should produce same result if ngenes, nevf and randseed are all the same)
#' @param modules a vector that sums up the ngenes, of the number of genes in each module
#' @param prob the probability that the effect size is not 0
#' @param sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from 
#' @return a list of 3 matrices, each of dimension ngenes * nevf
#' @examples 
#' GeneEffects()

GeneEffects_modules <- function(ngenes,nevf,modules,randseed,prob=1,sd){
	set.seed(randseed)	
	nmodules <- length(modules)
	if(nmodules >= nevf){stop(paste('not enough evfs to support',nmodules,'modules (has to at least have the same number,)'))}
	if(sum(modules) != ngenes){stop(paste('the number of genes in each module need to sum up to the number of total genes'))}
	set.seed(randseed)
	effect <- lapply(c('kon','koff','s'),function(param){
		effect1 <- lapply(c(1:ngenes),function(i){
			nonzero <- sample(size=nevf-nmodules,x=c(0,1),prob=c((1-prob),prob),replace=T)
			nonzero[nonzero!=0]=rnorm(sum(nonzero),mean=0,sd=sd)
			return(nonzero)
		})
		effect1 <- do.call(rbind,effect1)
		start_module <- c(0,cumsum(modules)[1:(nmodules-1)])
		effect2 <- lapply(c(1:nmodules),function(i){
			mod_effect <- rep(0,ngenes)
			mod_effect[(start_module[i]+1):(start_module[i]+modules[i])] = rnorm(modules[i],mean=0,sd=sd)
			return(mod_effect)
		})
		effect2 <- do.call(cbind,effect2)
		cbind(effect1,effect2)
	})
	modules <- lapply(c(1:nmodules),function(i){rep(i,modules[i])})
	modules <- do.call(c,modules)
	return(list(effect,modules))
}


#' Getting GeneEffects matrices 
#'
#' This function randomly generates the effect size of each evf on the dynamic expression parameters
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param randomseed (should produce same result if ngenes, nevf and randseed are all the same)
#' @param modules a vector that sums up the ngenes, of the number of genes in each module
#' @param prob the probability that the effect size is not 0
#' @param sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from 
#' @return a list of 3 matrices, each of dimension ngenes * nevf
#' @examples 
#' GeneEffects()
GeneEffects_corr <- function(ngenes,nevf,modules,randseed,prob,sd){
	set.seed(randseed)
	nmodules <- length(modules)
	effect <- lapply(c('kon','koff','s'),function(param){
		module_mean <- lapply(c(1:nmodules),function(i){
			nonzero <- sample(size=nevf,x=c(0,1),prob=c((1-prob),prob),replace=T)
			nonzero[nonzero!=0]=rnorm(sum(nonzero),mean=0,sd=sd)
			return(nonzero)
		})
		effect <- lapply(c(1:nmodules),function(i){
			sapply(c(1:nevf),function(j){rnorm(modules[i],mean=module_mean[[i]][j],sd=sd)})
		})
		return(do.call(rbind,effect))
	})
	modules <- lapply(c(1:nmodules),function(i){rep(i,modules[i])})
	modules <- do.call(c,modules)
	return(list(effect,modules))
}

#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range by using first a logistic function and then adding/dividing by constant to their correct range
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param param_dist the fitted parameter distribution to sample from 
#' @param bimod the bimodality constant
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params <- function(gene_effects,evf,param_dist,bimod){
	nparams <- length(param_dist[,1])
	params <- lapply(gene_effects,function(X){X %*% evf})
	params <- do.call(cbind,params)
	params <- apply(params,2,function(x){1/(1+exp(-x))})
	params[,1] = sort(param_dist[,1])[ceiling(params[,1]*nparams)]
	params[,2] = sort(param_dist[,2])[ceiling(params[,2]*nparams)]
	params[,3] = sort(param_dist[,3])[ceiling(params[,3]*nparams)]
	# params[,1] <- (params[,1]*7-1)
	# params[,2] <- params[,2]*7-2
	# params[,3] <- 10^(params[,3]*3)
	params[,1] <- log(base=10,params[,1])
	params[,2] <- log(base=10,params[,2])
	params[,1] <- 10^(params[,1] - (params[,1]+0.5)*bimod)
	params[,2] <- 10^(params[,2] - (params[,2]+1)*bimod)
	return(params)
}

#' Sample read counts from master equation
#'
#' This function takes ngenes distributions of transcript counts calculated by the master equation and sample from them
#' @param sim_dist master equation
#' @param cells the number of times to sample from each distribution
#' @return counts is a matrix of dimension ngenes * ncells
#' @examples 
#' SampleExprs()
SampleExprs <- function(sim_dist,ncells){
	counts=lapply(sim_dist,function(X){
		sample((c(1:length(X))-1),prob=X,size=ncells,replace=T)
	})
	counts=do.call(rbind,counts)
	return(counts)
}
#' Sample from true read counts to mimic capture/RT efficiency (only part of the transcripts are represented in the library)
#'
#' This function takes the true transcript counts and samples from it using a binomial distribution
#' @param sim_exprs the true matrix of transcript counts with ngenes rows and ncells columns
#' @param alpha the mean efficiency of capture
#' @param alpha_sd the standard deviation of the efficieny of capture for different cells
#' @return dropped_exprs The expression after adjusting for library size (alpha)
#' @examples 
#' TrueCounts2Dropped()
TrueCounts2Dropped <- function(sim_exprs,alpha,alpha_sd){
	alphas <- rnorm(length(sim_exprs[1,]),mean=alpha,sd=alpha_sd)
	alphas[alphas<0] <- 0	
	alphas[alphas>1] <- 1
	dropped_exprs <- lapply(c(1:length(alphas)),function(i){
		X <- sim_exprs[,i]
		sapply(X,function(Y){rbinom(n=1,size=Y,prob=alphas[i])})
	})
	dropped_exprs <- do.call(cbind,dropped_exprs)
	return(dropped_exprs)
}
#' Add Bias Terms and batch effects on the dropped count matrix
#'
#' randomly sample GC and length content of the simulated genes, and bias the exrpession according to it. Also add batch effect that is applied to all the the cells in the same batch, and a random epsilon error term to every number  
#' @param dropped_exprs dropped expression matrix
#' @param nbins number of bins when simulating quantile bias
#' @param randseed random seed when assigning simulated genes to GC and length bins
#' @param gcslope the slope of the GC bias line (the middle category is 1)
#' @param lenslope the slope of the length bias line (middle category is 1)
#' @param batch a vector of length ngenes that increase or decreases the expression of the same genes in the same batch by X fold
#' @param epsilon The standard deviation of the log error term that is multiplied to each gene 
#' @return biased_exprs The exrpession matrix 
#' @examples 
#' Dropped2Biased()
Dropped2Biased <- function(dropped_exprs,nbins,randseed,gcslope,lenslope,batch,epsilon){
	if(gcslope/2>=1 | lenslope/2>=1){
		print('ERROR:the two slope variables must be smaller than 2')
		return(NA)
	}
	gcslope <- gcslope/nbins; lenslope <- lenslope/nbins
	ngenes<- length(dropped_exprs[,1])
	set.seed(randseed)
	rand_gene1 <- sample(c(1:nbins),ngenes,replace=T)
	rand_gene2 <- sample(c(1:nbins),ngenes,replace=T)
	gc_bias <- c(1:nbins)*gcslope
	gc_bias<-(gc_bias+(1-median(gc_bias)))[rand_gene1]
#gc_bias and len_bias does not change gene-gene correlations because they are constant across cells
	len_bias <- (-c(1:nbins))*lenslope
	len_bias<-(len_bias+(1-median(len_bias)))[rand_gene2]
	biased_exprs <- apply(dropped_exprs,2,function(X){
		jiggle <- rnorm(length(X),mean=0,sd=epsilon)
		Y <- X * gc_bias * len_bias * batch * exp(jiggle)
		return(Y)
	})
	return(list(biased_exprs, list(gc_bias, len_bias, batch)))
}

#' Simulate 1 homogenous population in 1 batch
#' 
#' This function takes a gene effect matrix, the parameters for sampling extrinsic, intronsic and experimental variation, and outputs results of simulation
#' @param evf_mean a vector of length nevf, the means of the normal distribution to sample evf from 
#' @param evf_sd a vector of length nevf, the standard deviation of the normal distribution to sample evf from
#' @param ncells a integer of the number of cells to simulate
#' @param randseed the random seed to generate the gc and length bins (should keep constant for each experiment regardless of the batch and population, but should change for replicates of simulations)
#' @param gene_effects a matrix of ngenes * nevf, generated by GeneEffects function
#' @param bimod the proportion of distance to move a point in the kon-koff space towards the center of the kon-koff space where the gene expression is bimodal
#' @param alpha library prep capture efficiency, or dropout of transcript: the number of of captured transcript is distributed as the number of successes with probability of success alpha in X trials, where X is the true number of transcripts
#' @param alpha_sd the standard deviation of alpha (to add nsimoise for different efficiency for each cell)
#' @param nbins the number of bins for gc and length bias
#' @param gcbias the magnitude of gc bias (a number between 0 and 2) 
#' @param lenbias the magnitude of length bias (a number between 0 and 2)
#' @param batch_mean The mean batch effect from which the batch effect of each gene is sampled from
#' @param batch_sd The standard deviation of batch effect from which the batch effect of each gene is sampled from
#' @param noise the standard deviation of a normal distribution where the log(noise) is sampled from
#' @return list a list of 4 elements: evf, true counts, sampled counts and biased counts. 
sim1Pop1Batch <- function(evf_mean, evf_sd,ncells,randseed,gene_effects,
	bimod,alpha,alpha_sd,nbins,gcbias,lenbias,batch,noise,matched_params){
	a=as.numeric(Sys.time())
	set.seed(a)
	ngenes <- length(gene_effects[[1]][,1])
	if(gcbias>=2 | lenbias>=2){return('Error: gcbias or lenbias must be smaller than 2')}
	evfs <- lapply(c(1:ncells),function(celli){
		evf <- sapply(c(1:length(evf_mean)),function(evfi){rnorm(1,evf_mean[evfi],evf_sd[evfi])})
		return(evf)
	})
	true_counts <- lapply(c(1:ncells),function(celli){
		true_counts <- EVF2TrueCounts(allparams,matched_params,sim_master,evfs[[celli]],gene_effects,bimod)
		return(true_counts)
	})
	true_counts <- do.call(cbind,true_counts)
	sampled_counts <- TrueCounts2Dropped(true_counts,alpha,alpha_sd)
	biased_counts <- Dropped2Biased(sampled_counts,nbins,randseed,gcbias,lenbias,batch,noise)
	return(list(evfs,true_counts,sampled_counts,biased_counts[[1]],biased_counts[[2]]))
}
#' Simulate 1 homogenous population in 1 batch
#' 
#' This function takes a gene effect matrix, the parameters for sampling extrinsic, intronsic and experimental variation, and outputs results of simulation
#' @param evf_mean a vector of length nevf, the means of the normal distribution to sample evf from 
#' @param evf_sd a vector of length nevf, the standard deviation of the normal distribution to sample evf from
#' @param ncells a integer of the number of cells to simulate
#' @param randseed the random seed to generate the gc and length bins (should keep constant for each experiment regardless of the batch and population, but should change for replicates of simulations)
#' @param gene_effects a matrix of ngenes * nevf, generated by GeneEffects function
#' @param bimod the proportion of distance to move a point in the kon-koff space towards the center of the kon-koff space where the gene expression is bimodal
#' @param alpha library prep capture efficiency, or dropout of transcript: the number of of captured transcript is distributed as the number of successes with probability of success alpha in X trials, where X is the true number of transcripts
#' @param alpha_sd the standard deviation of alpha (to add nsimoise for different efficiency for each cell)
#' @param nbins the number of bins for gc and length bias
#' @param gcbias the magnitude of gc bias (a number between 0 and 2) 
#' @param lenbias the magnitude of length bias (a number between 0 and 2)
#' @param batch_mean The mean batch effect from which the batch effect of each gene is sampled from
#' @param batch_sd The standard deviation of batch effect from which the batch effect of each gene is sampled from
#' @param noise the standard deviation of a normal distribution where the log(noise) is sampled from
#' @return list a list of 4 elements: evf, true counts, sampled counts and biased counts. 
simCont <- function(evf_mean, evf_sd,ncells,randseed,gene_effects,
	bimod,alpha,alpha_sd,nbins,gcbias,lenbias,batch,noise,matched_params){
	a=as.numeric(Sys.time())
	set.seed(a)
	ngenes <- length(gene_effects[[1]][,1])
	if(gcbias>=2 | lenbias>=2){return('Error: gcbias or lenbias must be smaller than 2')}
	evfs <- lapply(c(1:ncells),function(celli){
		evf <- sapply(c(1:length(evf_mean)),function(evfi){rnorm(1,evf_mean[evfi],evf_sd[evfi])})
		return(evf)
	})
	true_counts <- lapply(c(1:ncells),function(celli){
		true_counts <- EVF2TrueCounts(allparams,matched_params,sim_master,evfs[[celli]],gene_effects,bimod)
		return(true_counts)
	})
	true_counts <- do.call(cbind,true_counts)
	sampled_counts <- TrueCounts2Dropped(true_counts,alpha,alpha_sd)
	biased_counts <- Dropped2Biased(sampled_counts,nbins,randseed,gcbias,lenbias,batch,noise)
	return(list(evfs,true_counts,sampled_counts,biased_counts[[1]],biased_counts[[2]]))
}
#' Simulate multiple discrete population in 1 batch
#' 
#' This function generates a evf variance-covariance matrix from a tree structure in the 
#' @param phyla a tree in the phylo format
#' @param nevf number of evfs 
#' @param nbatch number of batches to simulate
#' @param evf_sd a vector of length nevf, the standard deviation of the normal distribution to sample evf from
#' @param ncells a vector of length npop, the number of cells in each population
#' @param randseed the random seed to generate the gc and length bins (should keep constant for each experiment regardless of the batch and population, but should change for replicates of simulations)
#' @param gene_effects a matrix of ngenes * nevf, generated by GeneEffects function
#' @param bimod the proportion of distance to move a point in the kon-koff space towards the center of the kon-koff space where the gene expression is bimodal
#' @param alpha a vector of length nbatches, library prep capture efficiency, or dropout of transcript: the number of of captured transcript is distributed as the number of successes with probability of success alpha in X trials, where X is the true number of transcripts
#' @param alpha_sd a vector of length nbatches, the standard deviation of alpha (to add nsimoise for different efficiency for each cell)
#' @param nbins the number of bins for gc and length bias
#' @param gcbias a vector of length nbatches, the magnitude of gc bias (a number between 0 and 2) 
#' @param lenbias a vector of length nbatches, the magnitude of length bias (a number between 0 and 2)
#' @param batch_mean a vector of length nbatch that specify the magnitude of batch effect for each batch
#' @param noise the standard deviation of a normal distribution where the log(noise) is sampled from
#' @return list a list of 2 elements: all_counts and meta. all_counts is a list of 3 matrices of nrow=number of genes and ncol=sum of number of cells in each population. The fist matrix is the real counts, the second matrix is the counts after dropout, and the third matrix is the biased counts. each row of meta correspond to 1 cell in the count matrix, and the most important information in it is the 'pop' column that tells which population a cell comes from

Npop1Batch <- function(phyla, nevf,
		evf_sd, ncells, 
	randseed, gene_effects, bimod, alpha, 
	alpha_sd,nbins,gcbias,lenbias, batch, noise, pop_evf_mean=NA){

	npop <- length(phyla$tip.label)
	if(length(ncells)==1){ncells <- rep(ncells,npop)}
	if(npop!=length(ncells)){
		print('number of populations specified by cell size vector has to be the same as the number of tips in the tree')
		stop()
	}
	if(is.na(pop_evf_mean)[1]){
		cor_evf_mean<-vcv.phylo(phyla,cor=T)
		pop_evf_mean<-mvrnorm(nevf,rep(0,npop),cor_evf_mean)		
	}
	ngenes <- length(gene_effects[[1]][,1])
   	batch <- exp(rnorm(ngenes,mean=0,sd=batch))

	#need to specify the cell population name for cell size because the cor_evf_mean is scrambled 
	pop_change <- lapply(c(1:npop),function(pop){
		evf_mean <- pop_evf_mean[,pop]
		sim1Pop1Batch(evf_mean=evf_mean, evf_sd=rep(evf_sd,nevf),ncells=ncells[pop],
			randseed=randseed,gene_effects=gene_effects,bimod=bimod,alpha=alpha,
			alpha_sd=alpha_sd,nbins=nbins,gcbias=gcbias,lenbias=lenbias,
			batch=batch,noise=noise)
	})

	all_counts<-lapply(c(1:3),function(i){
		do.call(cbind,lapply(pop_change,function(X){X[[i+1]]}))	
	})
	bias <- lapply(pop_change,function(X){X[[5]]})
	evfs <- do.call(cbind,lapply(pop_change,function(X){X[[1]]}))
	meta <- data.frame(beta=rep(bimod,sum(ncells)),sigma=rep(evf_sd,sum(ncells)),alpha=rep(alpha,sum(ncells)),
	pop=do.call(c,lapply(c(1:length(ncells)),function(x){rep(colnames(cor_evf_mean)[x],ncells[x])})))
	return(list(all_counts,meta,evfs,bias))
}

#' Simulate multiple discrete population in multiple batches
#' 
#' This function generates a evf variance-covariance matrix from a tree structure in the 
#' @param phyla a tree in the phylo format
#' @param nevf number of evfs 
#' @param nbatches number of batches
#' @param evf_sd a vector of length nevf, the standard deviation of the normal distribution to sample evf from
#' @param ncells a matrix of size nbatch*npop, each row is a vector of length npop, and represents the number of cells from each population in each batch
#' @param randseed the random seed to generate the gc and length bins (should keep constant for each experiment regardless of the batch and population, but should change for replicates of simulations)
#' @param gene_effects a matrix of ngenes * nevf, generated by GeneEffects function
#' @param bimod the proportion of distance to move a point in the kon-koff space towards the center of the kon-koff space where the gene expression is bimodal
#' @param alpha a vector of length nbatches: library prep capture efficiency, or dropout of transcript: the number of of captured transcript is distributed as the number of successes with probability of success alpha in X trials, where X is the true number of transcripts
#' @param alpha_sd a vector of length nbatches: the standard deviation of alpha (to add nsimoise for different efficiency for each cell)
#' @param nbins the number of bins for gc and length bias
#' @param gcbias a vector of length nbatches: the magnitude of gc bias (a number between 0 and 2) 
#' @param lenbias a vector of length nbatches: the magnitude of length bias (a number between 0 and 2)
#' @param batch a vector of length nbatches: the constant to multiply all transcripts by 
#' @param noise a vector of length nbatches: the standard deviation of a normal distribution where the log(noise) is sampled from
#' @return list a list of 2 elements: all_counts and meta. all_counts is a list of 3 matrices of nrow=number of genes and ncol=sum of number of cells in each population. The fist matrix is the real counts, the second matrix is the counts after dropout, and the third matrix is the biased counts. each row of meta correspond to 1 cell in the count matrix, and the most important information in it is the 'pop' column that tells which population a cell comes from and 'batch' the column that tells which batch a cell is from

NpopNBatch <- function(phyla, nevf,nbatches,
	evf_sd, ncells, 
	randseed, gene_effects, bimod, alpha, 
	alpha_sd,nbins,gcbias,lenbias, batch, noise,
	same_evf_mean){
	if(!(length(ncells[,1])==nbatches & length(alpha)==nbatches &
		length(alpha_sd)==nbatches & length(gcbias)==nbatches &
		length(lenbias)==nbatches & length(batch)==nbatches & length(noise)==nbatches)){
		print('number of rows in ncells, length of alpha, alpha_sd, gcbias, lenbias, batch and noise has to be the same as the number of batches')
		stop()
	}
	npop <- length(phyla$tip.label)
	if(same_evf_mean==T){
		cor_evf_mean<-vcv.phylo(phyla,cor=T)
		pop_evf_mean<-mvrnorm(nevf,rep(0,npop),cor_evf_mean)		
	}else(pop_evf_mean=NA)
	data<-lapply(c(1:nbatches),function(i){
		result<-Npop1Batch(phyla=phyla,nevf=nevf,
			evf_sd=evf_sd, ncells=ncells[i,],
			randseed=randseed,gene_effects=gene_effects,bimod=bimod,alpha=alpha[i],
			alpha_sd=alpha_sd[i],nbins=nbins,
			gcbias=gcbias[i],lenbias=lenbias[i],batch=batch[i], noise=noise[i],pop_evf_mean=pop_evf_mean)
		result[[2]]$batch <- rep(i,sum(ncells[i,]))
		return(result)
	})
	all_counts <- lapply(c(1:3),function(i){
		do.call(cbind,lapply(data,function(X){X[[1]][[i]]}))
	})
	meta <- do.call(rbind,lapply(data,function(X){X[[2]]}))
	evfs <- lapply(data,function(X){X[[3]]})
	#batch, pop, evf
	bias <- lapply(data,function(X){X[[4]]})
	#batch, pop, gc, len, batch
	return(list(all_counts,meta,evfs,bias))
}


#' Test results of scaling normalization
#' 
#' This function runs rpm, deseq, tmm, uq (quantile normalization), scran, and clt normalization
#' @param sim_results output of sim1Pop1Batch, list a list of 4 elements: evf, true counts, sampled counts and biased counts.
#' @return cor_compare a matrix of ngenes * 2*(1+nmethods), each row is a gene, and the columns are the correlation between the true counts to dropped counts, dropped counts corrected, biased counts and dropped counts corrected


### Compute value of impulse model

#' Compute value of impulse function given parameters.
#' 
#' Compute value of impulse function given parameters.
#' Enforces lower bound on value of function to avoid numerical
#' errors during model fitting.
#' 
#' @seealso Compiled version: \link{evalImpulse_comp}
#' 
#' @param vecImpulseParam (numeric vector number of impulse model parameters)
#' \{beta, h0, h1, h2, t1, t2\}
#' Vector of impulse model parameters.
#' @param vecTimepoints (numeric vector length number of time points) 
#' Time points to be evaluated.
#' 
#' @return vecImpulseValue (vec number of vecTimepoints) 
#' Model values for given time points.
#' 
#' @author David Sebastian Fischer
evalImpulse <- function(vecImpulseParam, vecTimepoints) {
    # beta is vecImpulseParam[1] h0 is vecImpulseParam[2] h1 is
    # vecImpulseParam[3] h2 is vecImpulseParam[4] t1 is vecImpulseParam[5]
    # t2 is vecImpulseParam[6]
    vecImpulseValue <- sapply(vecTimepoints, function(t) {
        (1/vecImpulseParam[3]) * 
            (vecImpulseParam[2] + (vecImpulseParam[3] - vecImpulseParam[2]) *
                 (1/(1 + exp(-vecImpulseParam[1] * (t - vecImpulseParam[5]))))) *
            (vecImpulseParam[4] + (vecImpulseParam[3] - vecImpulseParam[4]) *
                 (1/(1 + exp(vecImpulseParam[1] * (t - vecImpulseParam[6])))))
    })
    vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)
    
    return(vecImpulseValue)
}

#' find closest matching simulated parameter sets 
#'
#' using kd-tree method find the simulated distribution that is closest to the observed distribution, and as diagnostic, plot the expression heat map of the true counts and the best-match simulated counts, and also the distribution of the parameters
#' @param sim_master the output of the master equation simulation serving as the database
#' @param allparams the parameters that produced the sim_master
#' @param log_count_bin the log scale bins
#' @param counts the real expression matrix that we are trying to fit to
#' @param plotting (boolean) whether to generate plots or not, default is true
#' @param samplename the prefix of the output plots 
#' @return match_params the best matched parameters
MatchParams <- function(sim_master,allparams,log_count_bin,counts,plotting=TRUE,samplename){
	simlogdist <- Sim_LogDist(sim_master,log_count_bin)
	truelogdist <- LogDist(counts, log_count_bin)
	outofrange <- is.na(rowMeans(truelogdist))
	counts <- counts[!outofrange,]
	truelogdist <- truelogdist[!outofrange,]
	#there are 6 genes ACTB (beta actin), GAPDH, Il17f(Interleukin 17f), Rpl41(Ribosomal Protein L41) , Il9(Interleukin 9), Rn18s.rs5(18s RNA) taht have greater than 10^4 expression and are ignore 
	best_matches=knnx.index(data=simlogdist,query=truelogdist,k=10,algorithm=c('kd_tree'))
	matchedlogdist=lapply(best_matches[,1],function(X){simlogdist[X,]})
	matchedlogdist=do.call(rbind,matchedlogdist)
	rownames(matchedlogdist)=rownames(truelogdist);colnames(matchedlogdist)=colnames(truelogdist)
	match_params <- allparams[best_matches[,1],] 
	if(plotting==T){
		p1 <- PlotCountHeatmap(truelogdist,rowMeans(counts),given_ord=NA,0.99,filename=paste(samplename,'.true.logged.jpeg',sep=''))
		p2 <- PlotCountHeatmap(matchedlogdist,rowMeans(counts),given_ord=NA,0.99,filename=paste(samplename,'.sim.logged.jpeg',sep=''))
		p3 <- PlotParamHist(match_params,samplename)
	}
	return(list(match_params,p1,p2,p3))
}

