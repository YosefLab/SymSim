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
EVF2TrueCounts <- function(allparams,sim_master,evf,gene_effects,bimod){
	nevf <- length(evf)
	ngenes <- length(gene_effects[[1]][,1])
	params <- Get_params(gene_effects,evf)
	params[,1] <- params[,1] - (params[,1]+0.5)*bimod
	params[,2] <- params[,2] - (params[,1]+1)*bimod
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
#' @param random seed (should produce same result if ngenes, nevf and randseed are all the same)
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
#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range by using first a logistic function and then adding/dividing by constant to their correct range
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params <- function(gene_effects,evf){
	params <- lapply(gene_effects,function(X){X %*% evf})
	params <- do.call(cbind,params)
	params <- apply(params,2,function(x){1/(1+exp(-x))})
	params[,1] <- 10^(params[,1]*7-1)
	params[,2] <- 10^(params[,2]*7-2)
	params[,3] <- 10^(params[,3]*3)
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
#' @param batch a constant that increase or decreases the expression in the same batch by X fold
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
	return(biased_exprs)
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
#' @param alpha_sd the standard deviation of alpha (to add noise for different efficiency for each cell)
#' @param nbins the number of bins for gc and length bias
#' @param gcbias the magnitude of gc bias (a number between 0 and 2) 
#' @param lenbias the magnitude of length bias (a number between 0 and 2)
#' @param batch the constant to multiply all transcripts by 
#' @param noise the standard deviation of a normal distribution where the log(noise) is sampled from
#' @return list a list of 4 elements: evf, true counts, sampled counts and biased counts. 
sim1Pop1Batch <- function(evf_mean, evf_sd,ncells,randseed,gene_effects,
	bimod,alpha,alpha_sd,nbins,gcbias,lenbias,batch,noise){
	if(gcbias>=2 | lenbias>=2){return('Error: gcbias or lenbias must be smaller than 2')}
	evf_mean <- rep(0,nevf)
	evf_sd <- rep(1,nevf)
	evfs <- lapply(c(1:ncells),function(celli){
		evf <- rnorm(evf_mean,evf_sd)
		return(evf)
	})
	true_counts <- lapply(c(1:ncells),function(celli){
		true_counts <- EVF2TrueCounts(allparams,sim_master,evfs[[celli]],gene_effects,bimod)
		return(true_counts)
	})
	true_counts <- do.call(cbind,true_counts)
	sampled_counts <- TrueCounts2Dropped(true_counts,alpha,alpha_sd)
	biased_counts <- Dropped2Biased(sampled_counts,nbins,randseed,gcbias,lenbias,batch,noise)
	return(list(evfs,true_counts,sampled_counts,biased_counts))
}
