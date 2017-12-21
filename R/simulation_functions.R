#############################################################
# Master Equation Related Functions
#############################################################

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
EVF2TrueCounts <- function(allparams,param_realdata,sim_master,evf,gene_effects,bimod, scale_s){
  if (missing(scale_s)) {scale_s <- 1}
	nevf <- length(evf) # this evf is for only one cell
	ngenes <- length(gene_effects[[1]][,1])
	params <- Get_params(gene_effects,evf,param_realdata,bimod, scale_s)
	best_matches <- FNN::knnx.index(data=allparams,query=params,k=1,algorithm=c('kd_tree'))
	sim <- sim_master[best_matches]
	sim_exprs <- SampleExprs(sim,1)
	return(sim_exprs) # sim_exprs is gene expression values in only one cell
}


#' Getting GeneEffects matrices 
#'
#' This function randomly generates the effect size of each evf on the dynamic expression parameters
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param randomseed (should produce same result if ngenes, nevf and randseed are all the same)
#' @param prob the probability that the effect size is not 0
#' @param geffect_mean the mean of the normal distribution where the non-zero effect sizes are dropped from 
#' @param geffect_sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from 
#' @return a list of 3 matrices, each of dimension ngenes * nevf
#' @examples 
#' GeneEffects()
GeneEffects <- function(ngenes,nevf,randseed,prob,geffect_mean,geffect_sd){
	set.seed(randseed)
	lapply(c('kon','koff','s'),function(param){
		effect <- lapply(c(1:ngenes),function(i){
			nonzero <- sample(size=nevf,x=c(0,1),prob=c((1-prob),prob),replace=T)
			nonzero[nonzero!=0]=rnorm(sum(nonzero),mean=geffect_mean,sd=geffect_sd)
			return(nonzero)
		})
		return(do.call(rbind,effect))
	})
}

#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range 
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), 
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param param_realdata the fitted parameter distribution to sample from 
#' @param bimod the bimodality constant
#' @param scale_s transcription rate will be multiplied by this factor to increase cell size
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
# param_realdata is the parameter distribution from real data
Get_params <- function(gene_effects,evf,param_realdata,bimod,scale_s){
  if (missing(scale_s)){scale_s <- 1}
  # generate parameters from evf and gene effects
  # for all genes in one cell
  params <- lapply(gene_effects,function(X){X %*% evf})
  params <- do.call(cbind,params)
  params <- apply(params,2,function(x){1/(1+exp(-x))}) 
  
  # map generated parameters to distribution of params from real data
  # values in params calculated from sigmoid function are between 0-1
  nparams <- length(param_realdata[,1])
  params[,1] = sort(param_realdata[,1])[ceiling(params[,1]*nparams)]
  params[,2] = sort(param_realdata[,2])[ceiling(params[,2]*nparams)]
  params[,3] = sort(param_realdata[,3])[ceiling(params[,3]*nparams)]
  
  # increase the bimodality of data by parameter beta
  params[,1] <- log(base=10,params[,1])
  params[,2] <- log(base=10,params[,2])
  params[,1] <- 10^(params[,1] - (params[,1]+0.5)*bimod)
  params[,2] <- 10^(params[,2] - (params[,2]+1)*bimod)
  params[,3] <- params[,3] * scale_s
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

#' This function simulates the amplification, library prep, and the sequencing processes.
#' @param true_counts_1cell the true transcript counts for one cell (one vector)
#' @param protocol a string, can be "ss2" or "10x"
#' @param rate_2cap the capture efficiency for this cell
#' @param gene_len gene lengths for the genes/transcripts, sampled from real human transcript length
#' @param amp_bias amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high
#' @param nPCR the number of PCR cycles
#' @param N_molecules_SEQ number of molecules sent for sequencing; sequencing depth
#' @return read counts (if protocol="ss2") or UMI counts (if protocol="10x)
#' @examples 
amplify_1cell <- function(true_counts_1cell, protocol, rate_2cap=0.1, gene_len, amp_bias, 
                          rate_2PCR=0.8, nPCR=18, N_molecules_SEQ){
  ngenes <- length(gene_len)
  if (protocol=="ss2"){load("SymSim/len2nfrag.RData")} else 
    if(protocol=="10x"){load("SymSim/len2prob3pri.RData")} # where should we keep the vairables?
  inds <- vector("list",2)
  # expand the original vector and apply capture efficiency
  # maintain a transcript index vector: which transcript the molecule belongs to
  expanded_res <- expand2binary(c(true_counts_1cell,1))
  expanded_vec <- expanded_res[[1]]; trans_idx <- expanded_res[[2]]
  
  inds[[1]] <- which(expanded_vec > 0); expanded_vec <- expanded_vec[inds[[1]]]
  trans_idx <- trans_idx[inds[[1]]]
  
  captured_vec <- expanded_vec; captured_vec[runif(length(captured_vec)) > rate_2cap] <- 0
  captured_vec[length(captured_vec)] <- 1
  inds[[2]] <- which(captured_vec > 0); captured_vec <- captured_vec[inds[[2]]]
  trans_idx <- trans_idx[inds[[2]]]
  
  amp_rate <- c((rate_2PCR+amp_bias[trans_idx[1:(length(trans_idx)-1)]]),1)
  temp <- runif(length(captured_vec)) < amp_rate
  temp <- temp*2+captured_vec-temp
  
  for (iPCR in 2:nPCR){
    eff <- runif(length(temp))*amp_rate
    v1 <- temp*(1-eff)
    round_down <- (v1-floor(v1)) < runif(length(v1))
    v1[round_down] <- floor(v1[round_down]); v1[!round_down] <- ceiling(v1[!round_down])
    temp <- v1 + 2*(temp-v1)
  }
  PCRed_vec <- temp
  
  if (protocol=="ss2"){ # add fragmentation step here
    temp_vec <- PCRed_vec
    for (i in seq(2,1,-1)){
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec; 
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec)-1)]
    amp_mol_count=numeric(ngenes);
    GI=c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell>0)){
      x=recovered_vec[(GI[i]+1):GI[i+1]]
      amp_mol_count[i]=sum(x)
    }
    
    # for every copy of each transcript, convert it into number of fragments
    frag_vec <- numeric(ngenes)
    for (igene in which(amp_mol_count>0)){
      frag_vec[igene] <- sum(sample(len2nfrag[as.character(gene_len[igene]),], 
                                    amp_mol_count[igene], replace = TRUE))}
    
    SEQ_efficiency=N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1) {read_count <- frag_vec} else{
      read_count <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)}) }
    return(read_count)
  } else if (protocol=="10x"){
    # fragmentation: 
    frag_vec <- sapply(1:(length(PCRed_vec)-1), function(igene)
    {return(rbinom(n=1, size = PCRed_vec[igene], 
                   prob = len2prob3pri[as.character(gene_len[trans_idx[igene]])] ))})
    SEQ_efficiency <- N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1){sequenced_vec <- frag_vec} else {
      sequenced_vec <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)})}
    
    temp_vec <- c(sequenced_vec,1)
    
    for (i in seq(2,1,-1)){
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec; 
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec)-1)]
    
    UMI=numeric(ngenes)
    GI=c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell>0)){
      x=recovered_vec[(GI[i]+1):GI[i+1]];
      UMI[i]=sum(x>0);
    }
    return(UMI)
  }
}



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
  print("function MatchParams() is used")
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


#' Creating an example tree with 5 tips
#' @param plotting True for plotting the tree on console, False for no plot 
#' @return a tree object
Phyla5 <- function(plotting=F){
	# par(mfrow=c(2,2))
	phyla<-rtree(2)
	phyla <- compute.brlen(phyla,1)
	tip<-rtree(2)
	tip <- compute.brlen(phyla,1)
	phyla<-bind.tree(phyla,tip,1)
	# if(plotting==T){plot(phyla)}
	phyla<-bind.tree(phyla,tip,2)
	# if(plotting==T){plot(phyla)}
	phyla<-bind.tree(phyla,tip,2)
	# if(plotting==T){plot(phyla)}
	phyla <- compute.brlen(phyla,c(1,1,1,1,1,0.2,0.2,3))
	edges <- cbind(phyla$edge,phyla$edge.length)
	edges <- cbind(c(1:length(edges[,1])),edges)
	connections <- table(c(edges[,2],edges[,3]))
	root <- as.numeric(names(connections)[connections==2])
	tips <- as.numeric(names(connections)[connections==1])
	phyla$tip.label <- as.character(tips)
	if(plotting==T){
	  plot(phyla,show.tip.label = F,lwd=2)
	  tiplabels(cex=2)
	  nodelabels(cex=2)
	}
	return(phyla)
}


#' Generating EVFs for cells sampled along the trajectory of cell development
#' @param phyla tree for cell developement
#' @param ncells number of cells
#' @param nevf1 Number of EVFs that do not have an impulse signal
#' @param nevf2 Number of EVFs with an impulse signal
#' @param tip The leaf that the path with impulse lead to
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree 
#' @param plotting Whether to plot the trajectory or not
#' @param plotname The 
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
ContinuousEVF <- function(phyla,ncells,nevf1,nevf2,tip,Sigma,plotting=T,plotname){
	edges <- cbind(phyla$edge,phyla$edge.length)
	edges <- cbind(c(1:length(edges[,1])),edges)
	connections <- table(c(edges[,2],edges[,3]))
	root <- as.numeric(names(connections)[connections==2])
	tips <- as.numeric(names(connections)[connections==1])
	internal <- as.numeric(names(connections)[connections==3])
	neutral <- SampleSubtree(root,0,0,edges,ncells,Sigma)
	neutral2 <- lapply(c(1:(nevf1-1)),function(evf_i){
		SampleSubtree(root,0,0,edges,ncells,Sigma,neutral=neutral)		
	})
	neutral2 <- lapply(neutral2,function(X){X[,4]})
	neutral2 <- do.call(cbind,neutral2)
 	pdf(file = plotname,width=15,height=5)
	evfs_wimpuls<- lapply(c(1:nevf2),function(evf_i){
	    impulse <-ImpulseEVFpertip(phyla, edges,root,tips,internal, neutral, tip,Sigma)
	    if(plotting==T){PlotRoot2Leave(impulse,tips,edges,root,internal)}
	    return(impulse)
  	})
  	dev.off()
  	re_order <- match(
  		apply(neutral[,c(1:3)],1,function(X){paste0(X,collapse='_')}),
  		apply(evfs_wimpuls[[1]][,c(1:3)],1,function(X){paste0(X,collapse='_')}))
  	evfs_wimpuls <- lapply(evfs_wimpuls,function(X){X[,4]})
  	evfs_wimpuls <- do.call(cbind,evfs_wimpuls)
  	evfs <- cbind(neutral[,4],neutral2,evfs_wimpuls)
  	colnames(evfs)<-c(rep('neutral',nevf1),rep('impulse',nevf2))
  	meta <- data.frame(pop=apply(neutral[,c(1:2)],1,function(X){paste0(X,collapse='_')}),depth=neutral[,3])
  	return(list(evfs[c(1:ncells),],meta[c(1:ncells),]))
}

#' Generating EVFs for cells sampled from tip populations from a tree
#' @param phyla tree for cell developement
#' @param ncells_total number of cells from all populations
#' @param min_popsize size of the rarest population
#' @param nevf Number of EVFs per cell
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree 
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the population each cell comes from (pop)
DiscreteEVF <- function(phyla, ncells_total, min_popsize, Sigma, nevf){
	npop <- length(phyla$tip.label) # number of populations
	
	# set the number of cells in each population
	# first give each population min_popsize cells
	# then randomly distribute the rest of cells to all populations except the smallest (first) one
	ncells_pop <- rep(min_popsize, npop)
	if (ncells_total <= min_popsize*npop) {
	  stop("The size of the smallest population is too big for the total number of cells")}
	temp <- sample(2:npop, (ncells_total-min_popsize*npop), replace = T)
	ncells_pop[2:npop] <- ncells_pop[2:npop] + table(temp)
	cor_evf_mean<-vcv.phylo(phyla,cor=T)
	pop_evf_mean<-mvrnorm(nevf,rep(0,npop),cor_evf_mean)	
	evfs <- lapply(c(1:npop),function(ipop){
		evf <- sapply(c(1:nevf),function(ievf){rnorm(ncells_pop[ipop],pop_evf_mean[ievf,ipop],Sigma)})
		return(evf)
	})
	evfs <- do.call(rbind,evfs)
	meta <- data.frame(pop=do.call(c,lapply(c(1:npop),function(i){rep(i,ncells_pop[i])})))
	return(list(evfs,meta))
}

#' Generate alpha vector for ncells
#' @param ncells number of cells
#' @param alpha the mean value of a normal distribution where alpha is sampled from
#' @param alpha_sd the standard deviation value of a normal distribution where alpha is sampled from
#' @return a vector with an alpha value for each cell
SampleAlpha <- function(ncells,alpha,alpha_sd){
	alphas <- rnorm(ncells,mean=alpha,sd=alpha_sd)
	alphas[alphas<0] <- 0	
	alphas[alphas>1] <- 1
	return(alphas)
}


#' Generate gene bias vector for each gene
#' @param ngenes number of gene 
#' @param nbatch number of batches 
#' @param batch_sd the standard deviation for the random per genebias in each batch (same for each batch)
#' @return a vector with an batch bias for each gene
SampleBatch <- function(ngenes,nbatch,batch_sd){
	batch <- lapply(c(1:nbatch),function(i){
	   	exp(rnorm(ngenes,mean=0,sd=batch_sd))
	})
	batch <- do.call(cbind,batch)
	return(batch)
}

#' Simulate count matrix from evf, gene effects and technical biases
#' @param ncells number of cells
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param nbatch number of batches
#' @param evfs matrix with ncells rows and nevf columns, can be generated by ContinuousEVF and DiscreteEVF
#' @param gene_effects a list of length 3, each element is a sparse matrix determining how each gene is affected by the EVF of that cell
#' @param alphas a vector of length ncells
#' @param beta a single float value between 0 and 1 determining the bimodality of entire dataset
#' @param gcbias a vector of length ngenes
#' @param lenbias a vector of length ngenes
#' @param batch a matrix of ngenes rows and nbatch columns
#' @param batch_id a vector of length ncells indicating which batch each cell comes from (randomly generated)
#' @param meta1 meta data generated from EVF function (population, depth in the tree, etc.)
#' @param epsilon a single float value the magnitude of independent noise for each gene in each cell
#' @return a list of 3 objects, one is a list of the count matrices (true counts, dropped counts and biased counts), the per cell meta data, and the per gene meta data
SimulateCounts <- function(ncells,ngenes,nevf,nbatch,evfs,gene_effects,
	alphas,beta,gcbias,lenbias,batch,batch_id,meta1,epsilon){
	#check dimension of input
	if(length(evfs[,1])!=ncells){print('ERROR: number of rows of evf is the number of cells');stop()}
	if(length(evfs[1,])!=nevf){print('ERROR: number of columns of evf is the number of evfs');stop()}
	if(length(alphas)!=ncells){print('ERROR:need one alpha for each cell');stop()}
	if(length(beta)!=1){print('ERROR:need only one beta for an experiment');stop()}
	if(length(gcbias)!=ngenes){print('ERROR: need GC bias for each gene');stop()}
	if(length(lenbias)!=ngenes){print('ERROR: need length bias for each gene');stop()}
	if(length(batch[,1])!=ngenes){print('ERROR: number of rows in batch effect is equal to the number of genes');stop()}
	if(length(batch[1,])!=nbatch){print('ERROR: number of columns if batch effect is equal to the number of bathces');stop()}
	
	
	biased_exprs <- lapply(c(1:ncells),function(i){
		X <- dropped_exprs[,i]
		jiggle <- rnorm(ngenes,mean=0,sd=epsilon)
		Y <- X * gcbias * lenbias * batch[,batch_id[i]] * exp(jiggle)
		return(Y)
	})
	
	biased_exprs <- do.call(cbind,biased_exprs)
	meta <- data.frame(batch=batch_id,alpha=alphas,evfs=evfs)
	meta <- cbind(meta1,meta)
	return(list(list(true_counts,dropped_exprs,biased_exprs),meta))
}

#' Simulate count matrix from evf, gene effects and technical biases
#' @param ncells number of cells
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param nbatch number of batches
#' @param evfs matrix with ncells rows and nevf columns, can be generated by ContinuousEVF and DiscreteEVF
#' @param gene_effects a list of length 3, each element is a sparse matrix determining how each gene is affected by the EVF of that cell
#' @param alphas a vector of length ncells
#' @param beta a single float value between 0 and 1 determining the bimodality of entire dataset
#' @param gcbias a vector of length ngenes
#' @param lenbias a vector of length ngenes
#' @param batch a matrix of ngenes rows and nbatch columns
#' @param batch_id a vector of length ncells indicating which batch each cell comes from (randomly generated)
#' @param epsilon a single float value the magnitude of independent noise for each gene in each cell

# this function needs substantial work!
# also need to make the size of clusters different
SimulateTrueCounts <- function(ncells_total,min_popsize=ncells_total,ngenes,
                            nevf=10,evf_type="one.population",Sigma=0.3,phyla=NULL,
                            gene_effects_sd=1,gene_effect_prob=0.3,
                            bimod=0.3,randseed=0){
  if(evf_type=='one.population'){
    evf_mean=rep(0,nevf); evf_sd=rep(Sigma,nevf)
    evfs <- lapply(c(1:ncells_total),function(celli){
      evf <- sapply(c(1:nevf),function(ievf){rnorm(1,evf_mean[ievf],evf_sd[ievf])})
      return(evf)
    })
    evf_res <- list(evfs=do.call(rbind, evfs), meta=data.frame(pop=rep(1, ncells_total)))
  } else if(evf_type=='discrete'){
    evf_res <- DiscreteEVF(phyla,ncells_total,min_popsize,Sigma,nevf)
  }else if(evf_type=='continuous'){
    evf_res <- ContinuousEVF(phyla,ncells_total,nevf1=nevf/2,nevf2=nevf/2,
                             tip=1,Sigma,plotting=T,plotname)		
  }
  set.seed(randseed)
  seed <- sample(c(1:1e5),size=3)
  gene_effects <- GeneEffects(ngenes=ngenes,nevf=nevf,randseed=seed[3],prob=gene_effect_prob, 
                              geffect_mean=0, geffect_sd=gene_effects_sd)
  true_counts <- lapply(c(1:ncells_total),function(celli){
    true_counts <- EVF2TrueCounts(allparams,param_realdata,sim_master,evf=evf_res[[1]][celli,],gene_effects,bimod=bimod, scale_s=1)
    return(true_counts)
  })
  true_counts <- do.call(cbind,true_counts)
  return(list(counts=true_counts, meta_cell=evf_res[[2]]))
}



True2ObservedCounts <- function(true_counts,meta_cell,nbatch=1,protocol,alpha_mean=0.1,alpha_sd=0.02,
                                lenslope=0.01,nbins=20,gene_len,amp_bias_limit=c(-0.2, 0.2),
                                batch_sd=0.5,epsilon=0.05,rate_2PCR=0.8,nPCR=16,
                                depth_mean, depth_sd,randseed=0){  
  ngenes <- dim(true_counts)[1]; ncells <- dim(true_counts)[2]
  amp_bias <- cal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
  
  rate_2cap_vec <- rnorm(ncells, mean = alpha_mean, sd=alpha_sd)
  rate_2cap_vec[which(rate_2cap_vec < 0.01)] <- 0.01
  depth_vec <- rnorm(ncells, mean = depth_mean, sd=depth_sd)
  depth_vec[which(depth_vec < 500)] <- 500
  
  observed_counts <- matrix(0, ngenes, ncells)
  for (icell in 1:ncells){
    observed_counts[, icell] <- amplify_1cell(true_counts_1cell =  true_counts[, icell], protocol=protocol, 
                                              rate_2cap=rate_2cap_vec[icell], gene_len=gene_len, amp_bias = amp_bias, 
                                              rate_2PCR=rate_2PCR, nPCR=nPCR, N_molecules_SEQ = depth_vec[icell]) 
  }
  meta_cell2 <- data.frame(alpha=rate_2cap_vec,depth=depth_vec)
  meta_cell <- cbind(meta_cell, meta_cell2)
  return(list(observed_counts, meta_cell))
}



