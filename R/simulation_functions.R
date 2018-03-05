
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

<<<<<<< 176a52e550dca83672076caa19801f69fc1b3d4b
=======
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

>>>>>>> update vignette files

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
#' sample from smoothed density function
#' @param nsample number of samples needed
#' @param den_fun density function estimated from density() from R default
SampleDen <- function(nsample,den_fun){
  probs <- den_fun$y/sum(den_fun$y)
  bw <- den_fun$x[2]-den_fun$x[1]
  bin_id <- sample(size=nsample,x=c(1:length(probs)),prob=probs,replace=T)
  counts <- table(bin_id)
  sampled_bins <- as.numeric(names(counts))
  samples <- lapply(c(1:length(counts)),function(j){
    runif(n=counts[j],min=(den_fun$x[sampled_bins[j]]-0.5*bw),max=(den_fun$x[sampled_bins[j]]+0.5*bw))
  })
  samples <- do.call(c,samples)
  return(samples)
}

#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range 
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), 
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param match_param_den the fitted parameter distribution density to sample from 
#' @param bimod the bimodality constant
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
<<<<<<< HEAD
<<<<<<< 176a52e550dca83672076caa19801f69fc1b3d4b
Get_params <- function(gene_effects,evf,match_params,bimod){
  nparams <- length(match_params[,1])
  params <- lapply(gene_effects,function(X){evf %*% t(X)})
=======
Get_params <- function(gene_effects,evf,match_param_den,bimod){
  params <- lapply(1:3, function(iparam){evf[[iparam]] %*% t(gene_effects[[iparam]])})
>>>>>>> de68f0d77d0fb7c5402e989a978af800e9338ce1
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    # X=matrix(data=c(1:10),ncol=2) 
    # this line is to check that the row and columns did not flip
    temp <- alply(X, 1, function(Y){Y})
    values <- do.call(c,temp)
    ranks <- rank(values)
    sorted <- sort(SampleDen(nsample=max(ranks),den_fun=match_param_den[[i]]))
    temp3 <- matrix(data=sorted[ranks],ncol=length(X[1,]),byrow=T)
    return(temp3)
  })
  
  bimod_perc <- 1
  ngenes <- dim(scaled_params[[1]])[2]; bimod_vec <- numeric(ngenes)
  bimod_vec[1:ceiling(ngenes*bimod_perc)] <- bimod
  bimod_vec <- c(rep(bimod, ngenes/2), rep(0, ngenes/2))
  scaled_params[[1]] <- apply(t(scaled_params[[1]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[2]] <- apply(t(scaled_params[[2]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[3]] <- t(apply(scaled_params[[3]],2,function(x){x<-abs(x)}))
  
  return(scaled_params)
}
#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range 
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), 
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param match_param_den the fitted parameter distribution density to sample from 
#' @param bimod the bimodality constant
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params3 <- function(gene_effects,evf,match_params,bimod){
  X <- evf[[1]] %*% t(gene_effects[[1]])
  temp <- alply(X, 1, function(Y){Y})
  values <- do.call(c,temp)
  ranks <- rank(values)
  PC1 <- prcomp(t(match_params))$rotation[,1]
  PC1 <- order(PC1)
  temp <- sample(c(1:length(PC1)),length(values),replace=T)
  temp <- sort(temp)
  PC1 <- PC1[temp]
  params <- match_params[PC1,]
  scaled_params <- lapply(c(1:3),function(i){
    matrix(data=params[ranks,i],ncol=length(X[1,]),byrow=T)
  })
  scaled_params[[1]]<-apply(scaled_params[[1]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[2]]<-apply(scaled_params[[2]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[3]]<-apply(scaled_params[[3]],2,function(x){x<-abs(x)})
  scaled_params <- lapply(scaled_params,t)
  return(scaled_params)
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

Get_params2 <- function(gene_effects,evf,bimod,ranges){
  params <- lapply(gene_effects,function(X){evf %*% t(X)})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    temp <- apply(X,2,function(x){1/(1+exp(-x))})
    temp2 <- temp*(ranges[[i]][2]-ranges[[i]][1])+ranges[[i]][1]
    return(temp2)
  })
  scaled_params[[1]]<-apply(scaled_params[[1]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[2]]<-apply(scaled_params[[2]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[3]]<-apply(scaled_params[[3]],2,function(x){x<-abs(x)})
  scaled_params <- lapply(scaled_params,t)
  return(scaled_params)
=======
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
>>>>>>> update vignette files
}


#' This function simulates the amplification, library prep, and the sequencing processes.
#' @param true_counts_1cell the true transcript counts for one cell (one vector)
#' @param protocol a string, can be "ss2" or "umi"
#' @param rate_2cap the capture efficiency for this cell
#' @param gene_len gene lengths for the genes/transcripts, sampled from real human transcript length
#' @param amp_bias amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high
#' @param nPCR the number of PCR cycles
#' @param N_molecules_SEQ number of molecules sent for sequencing; sequencing depth
#' @return read counts (if protocol="ss2") or UMI counts (if protocol="umi)
amplify_1cell <- function(true_counts_1cell, protocol, rate_2cap=0.1, gene_len, amp_bias, 
                          rate_2PCR=0.8, nPCR=18, N_molecules_SEQ){
  ngenes <- length(gene_len)
  if (protocol=="ss2"){load("SymSim/len2nfrag.RData")} else 
    if(protocol=="umi"){load("SymSim/len2prob3pri.RData")} # where should we keep the vairables?
  inds <- vector("list",2)
  # expand the original vector and apply capture efficiency
  # maintain a transcript index vector: which transcript the molecule belongs to
  expanded_res <- expand2binary(c(true_counts_1cell,1))
  expanded_vec <- expanded_res[[1]]; trans_idx <- expanded_res[[2]]
  
  inds[[1]] <- which(expanded_vec > 0); expanded_vec <- expanded_vec[inds[[1]]]
  trans_idx <- trans_idx[inds[[1]]]
  
  captured_vec <- expanded_vec; captured_vec[runif(length(captured_vec)) > rate_2cap] <- 0
  if (sum(captured_vec[1:(length(captured_vec)-1)]) < 1) {return(rep(0, ngenes))}
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
    if (any(is.na(round_down))) {stop(sprintf("there is NA values in round_down, capture rate=%4.2f, total true counts=%d", 
                                              rate_2cap, sum(true_counts_1cell)))}
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
  } else if (protocol=="umi"){
    # fragmentation: 
    frag_vec <- sapply(1:(length(PCRed_vec)-1), function(igene)
    {return(rbinom(n=1, size = PCRed_vec[igene], 
                   prob = len2prob3pri[as.character(gene_len[trans_idx[igene]])] ))})
    # another 10 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:3){
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
    }
    for (iPCR in 4:8){
      frag_vec <- frag_vec + round(frag_vec*rate_2PCR)
    }
    
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

<<<<<<< 176a52e550dca83672076caa19801f69fc1b3d4b
=======
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
>>>>>>> update vignette files


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

Phyla3 <- function(plotting=F){
  # par(mfrow=c(2,2))
  phyla<-rtree(2)
  phyla <- compute.brlen(phyla,1)
  tip<-rtree(2)
  tip <- compute.brlen(phyla,1)
  phyla<-bind.tree(phyla,tip,1)
  phyla <- compute.brlen(phyla,c(1,1,1,2))
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
#' @param n_nd_evf Number of EVFs that do not have an impulse signal
#' @param n_de_evf Number of EVFs with an impulse signal
#' @param tip The leaf that the path with impulse lead to
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree 
#' @param plotting Whether to plot the trajectory or not
#' @param plotname The 
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
ContinuousEVF <- function(phyla,ncells,n_nd_evf,n_de_evf,impulse=T,evf_center=1,vary='all',
                          Sigma,plotting=T,plotname='cont_evf.pdf',seed){
  set.seed(seed)
  edges <- cbind(phyla$edge,phyla$edge.length)
  edges <- cbind(c(1:length(edges[,1])),edges)
  edges[,4] <- edges[,4]/mean(vcv.phylo(phyla))
  connections <- table(c(edges[,2],edges[,3]))
  root <- as.numeric(names(connections)[connections==2])
  tips <- as.numeric(names(connections)[connections==1])
  internal <- as.numeric(names(connections)[connections==3])
  if(vary=='all'){
    N_DE_evfs =c(n_de_evf,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='kon'){
    N_DE_evfs =c(n_de_evf,0,0)    
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='koff'){
    N_DE_evfs =c(0,n_de_evf,0)    
    N_ND_evfs =c(n_de_evf+n_nd_evf,n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='s'){
    N_DE_evfs =c(0,0,n_de_evf)    
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf+n_de_evf,n_nd_evf)
  }else if(vary=='except_kon'){
    N_DE_evfs =c(0,n_de_evf,n_de_evf)    
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='except_koff'){
    N_DE_evfs =c(n_de_evf,0,n_de_evf)    
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_nd_evf)
  }else if(vary=='except_s'){
    N_DE_evfs =c(n_de_evf,n_de_evf,0)    
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf+n_de_evf)
  }
  neutral <- SampleSubtree(root,0,evf_center,edges,ncells,Sigma)
  param_names <- c("kon", "koff", "s")
  evfs <- lapply(c(1:3),function(parami){
    nd_evf <- lapply(c(1:N_ND_evfs[parami]),function(ievf){
      rnorm(ncells,evf_center,Sigma)
    })
    nd_evf <- do.call(cbind,nd_evf)
    if(N_DE_evfs[parami]!=0){
      #if there is more than 1 de_evfs for the parameter we are looking at
      if(impulse==T){
        pdf(file = plotname,width=15,height=5)
        tip <- rep(tips,ceiling(N_DE_evfs[parami]/length(tips)))
        de_evf <- lapply(c(1:N_DE_evfs[parami]),function(evf_i){
          impulse <-ImpulseEVFpertip(phyla, edges,root,tips,internal, neutral, tip[evf_i],Sigma,evf_center=evf_center)
          if(plotting==T){PlotRoot2Leave(impulse,tips,edges,root,internal)}
          re_order <- match(
            apply(neutral[,c(1:3)],1,function(X){paste0(X,collapse='_')}),
            apply(impulse[,c(1:3)],1,function(X){paste0(X,collapse='_')}))            
          return(impulse[re_order,])
        })
        dev.off()
      }else{
        de_evf <- lapply(c(1:N_DE_evfs[parami]),function(evf_i){
          SampleSubtree(root,0,evf_center,edges,ncells,Sigma*50,neutral=neutral)    
        })
      }
      de_evf <- lapply(de_evf,function(X){X[,4]})
      de_evf <- do.call(cbind,de_evf)
      de_evf <- de_evf[c(1:ncells),]
      evfs <- cbind(nd_evf,de_evf)
      colnames(evfs)<-c(
        paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_'),
        paste(param_names[parami],rep('DE',length(de_evf[1,])),c(1:length(de_evf[1,])),sep='_'))
    }else{
      evfs <- nd_evf
      colnames(evfs)<-paste(param_names[parami],rep('nonDE',length(nd_evf[1,])),c(1:length(nd_evf[1,])),sep='_')
    }
    return(evfs)
  })
  meta <- data.frame(pop=apply(neutral[,c(1:2)],1,function(X){paste0(X,collapse='_')}),depth=neutral[,3])
  return(list(evfs,meta[c(1:ncells),]))
  # note previously the number of sampled evfs and meta isn't necessarily ncells? 
}
#' Generating EVFs for cells sampled from tip populations from a tree
#' @param phyla tree for cell developement
#' @param ncells_total number of cells from all populations
#' @param min_popsize size of the rarest population
#' @param i_minpop to specify which population has the smallest size
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree 
#' @param nevf Number of EVFs per cell
#' @param evf_center the value used to generated evf means. Suggested value is 1
#' @param percent_DEevf the percentage of differential EVFs out of all evfs (nevf)
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the population each cell comes from (pop)
<<<<<<< HEAD
<<<<<<< 176a52e550dca83672076caa19801f69fc1b3d4b
DiscreteEVF <- function(phyla, ncells_total, min_popsize, Sigma, nevf,seed){
  set.seed(seed)
	npop <- length(phyla$tip.label) # number of populations
=======
DiscreteEVF <- function(phyla, ncells_total, min_popsize, Sigma, nevf){
	npop <- length(phyla$tip.label) # number of populations
	
>>>>>>> update vignette files
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
=======
DiscreteEVF <- function(phyla, ncells_total, min_popsize, i_minpop, Sigma, n_nd_evf, n_de_evf, vary, evf_center, seed){
  set.seed(seed)
  npop <- length(phyla$tip.label)
  # set the number of cells in each population: first give each population min_popsize cells
  # then randomly distribute the rest of cells to all populations except the smallest one
  ncells_pop <- rep(min_popsize, npop)
  if (ncells_total <= min_popsize*npop) {
    stop("The size of the smallest population is too big for the total number of cells")}
  larger_pops <- setdiff(1:npop, i_minpop)
  temp <- sample(larger_pops, (ncells_total-min_popsize*npop), replace = T)
  ncells_pop[larger_pops] <- ncells_pop[larger_pops] + table(temp)
  
  vcv_evf_mean <- vcv.phylo(phyla,cor=T)
  param_names <- c("kon", "koff", "s")
  if(vary=='all'){
    N_DE_evfs =c(n_de_evf,n_de_evf,n_de_evf)
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='kon'){
    N_DE_evfs =c(n_de_evf,0,0)    
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='koff'){
    N_DE_evfs =c(0,n_de_evf,0)    
    N_ND_evfs =c(n_de_evf+n_nd_evf,n_nd_evf,n_de_evf+n_nd_evf)
  }else if(vary=='s'){
    N_DE_evfs =c(0,0,n_de_evf)    
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf+n_de_evf,n_nd_evf)
  }else if(vary=='except_kon'){
    N_DE_evfs =c(0,n_de_evf,n_de_evf)    
    N_ND_evfs =c(n_nd_evf+n_de_evf,n_nd_evf,n_nd_evf)
  }else if(vary=='except_koff'){
    N_DE_evfs =c(n_de_evf,0,n_de_evf)    
    N_ND_evfs =c(n_nd_evf,n_de_evf+n_nd_evf,n_nd_evf)
  }else if(vary=='except_s'){
    N_DE_evfs =c(n_de_evf,n_de_evf,0)    
    N_ND_evfs =c(n_nd_evf,n_nd_evf,n_nd_evf+n_de_evf)
  }
  
  if (sum(N_DE_evfs) < 5) {warning("The number of DE evfs is less than 5; in the case of a small number of DE evfs, the structure of generated data 
	                       may not closely follow the input tree. One can either increase nevf or percent_DEevf to avoid this warning.")}
  
  evfs <- lapply(1:3, function(iparam){
    if (N_ND_evfs[iparam] > 0) {
      pop_evf_nonDE <- lapply(c(1:npop),function(ipop){
        evf <- sapply(c(1:(N_ND_evfs[iparam])),function(ievf){rnorm(ncells_pop[ipop],evf_center,Sigma)})
        return(evf)
      })
      pop_evf_nonDE <- do.call(rbind,pop_evf_nonDE)
      colnames(pop_evf_nonDE) <- rep('nonDE',N_ND_evfs[iparam])
    } else {pop_evf_nonDE <- NULL}
    if (N_DE_evfs[iparam] > 0){
      pop_evf_mean_DE <- mvrnorm(N_DE_evfs[iparam],rep(evf_center,npop),vcv_evf_mean)
      pop_evf_DE <- lapply(c(1:npop),function(ipop){
        evf <- sapply(c(1:N_DE_evfs[iparam]),function(ievf){rnorm(ncells_pop[ipop],pop_evf_mean_DE[ievf,ipop],Sigma)})
        return(evf)
      })
      pop_evf_DE <- do.call(rbind,pop_evf_DE)
      colnames(pop_evf_DE) <- rep('DE',N_DE_evfs[iparam])
    } else {pop_evf_DE <- NULL}
    
    evfs_per_param <- cbind(pop_evf_nonDE,pop_evf_DE)
    colnames(evfs_per_param) <- sprintf("%s_%s_evf%d", param_names[iparam],colnames(evfs_per_param), 
                                        1:(N_ND_evfs[iparam]+N_DE_evfs[iparam]))
    return(evfs_per_param)
  })
  meta <- data.frame(pop=do.call(c,lapply(c(1:npop),function(i){rep(i,ncells_pop[i])})))
  return(list(evfs,meta))
>>>>>>> de68f0d77d0fb7c5402e989a978af800e9338ce1
}
<<<<<<< HEAD


<<<<<<< HEAD
<<<<<<< 176a52e550dca83672076caa19801f69fc1b3d4b
=======
=======
>>>>>>> de68f0d77d0fb7c5402e989a978af800e9338ce1

=======
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
>>>>>>> update vignette files

>>>>>>> b1c727f6ac83412053a415eccd221e3a5179bbd7
#' Generate both evf and gene effect and simulate true transcript counts
#' @param ncells_total number of cells
#' @param min_popsize the number of cells 
#' @param i_minpop specifies which population has the smallest size
#' @param ngenes number of genes
#' @param evf_center the value which evf mean is generated from
#' @param nevf number of evfs
<<<<<<< 176a52e550dca83672076caa19801f69fc1b3d4b
#' @param evf_type string that is one of the following: 'one.population','discrete','continuous'
#' @param percent_DEevf percentage of differential evfs between populations
#' @param Sigma parameter of the std of evf values within the same population
#' @param phyla the cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly or read from newick format file using the ape package
#' @param gene_effect_prob the probability that the effect size is not 0
#' @param gene_effect_sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from 
#' @param match_params_den empirical density function of the kon,koff and s parameter estimated from real data
#' @param bimod the amount of increased bimodality in the transcript distribution, 0 being not changed from the results calculated using evf and gene effects, and 1 being all genes are bimodal
#' @param randseed random seed
#' @param SE return summerized experiment rather than a list of elements, default is False
#' @return a list of 4 elements, the first element is true counts, second is the gene level meta information, the third is cell level meta information, including a matrix of evf and a vector of cell identity, and the fourth is the parameters kon, koff and s used to simulation the true counts
SimulateTrueCounts <- function(ncells_total,min_popsize,i_minpop=1,ngenes, 
                               evf_center=1,evf_type="one.population",nevf=10,
<<<<<<< HEAD
                               n_de_evf=0,impulse=T,vary='all',
                               Sigma=0.5,phyla=NULL,geffect_mean=0,gene_effects_sd=1,gene_effect_prob=0.3,
                               bimod=0.2,param_realdata="zeisel.imputed",joint=F,randseed,SE=F){
  set.seed(randseed)
  n_nd_evf=nevf-n_de_evf
=======
                               n_nd_evf=nevf,n_de_evf=0,impulse=T,vary='all',
                               Sigma=0.5,phyla=NULL,geffect_mean=0,gene_effects_sd=1,gene_effect_prob=0.3,
                               bimod=0.2,param_realdata="zeisel.imputed",joint=F,randseed,SE=F){
  set.seed(randseed)
>>>>>>> b1c727f6ac83412053a415eccd221e3a5179bbd7
  seed <- sample(c(1:1e5),size=2)
  param_names <- c("kon", "koff", "s")
  if(evf_type=='one.population'){
    evf_mean=rep(evf_center,nevf); evf_sd=rep(Sigma,nevf)
    evfs <- lapply(1:3, function(iparam){
      evfs_per_param <- lapply(c(1:ncells_total),function(celli){
        evf <- sapply(c(1:nevf),function(ievf){rnorm(1,evf_mean[ievf],evf_sd[ievf])})
        return(evf)})
      evfs_per_param <- do.call(rbind, evfs_per_param)
      colnames(evfs_per_param) <- sprintf("%s_evf%d", param_names[iparam], 1:nevf)
      return(evfs_per_param)})
    evf_res <- list(evfs=evfs, meta=data.frame(pop=rep(1, ncells_total)))
  } else if(evf_type=='discrete'){
    evf_res <- DiscreteEVF(phyla,ncells_total,min_popsize,i_minpop=i_minpop,Sigma,
                           n_nd_evf, n_de_evf, vary=vary, evf_center=evf_center, seed=seed[1])
  }else if(evf_type=='continuous'){
    n_de <- round(nevf*percent_DEevf)
    evf_res <- ContinuousEVF(phyla,ncells_total,n_nd_evf=nevf-n_de,n_de_evf=n_de,
                             evf_center=evf_center,vary=vary,impulse=impulse,
                             Sigma,plotting=T,seed=seed[1])    
  }
  gene_effects <- GeneEffects(ngenes=ngenes,nevf=nevf,randseed=seed[2],prob=gene_effect_prob,
                              geffect_mean=geffect_mean,geffect_sd=gene_effects_sd)
  if(!is.null(param_realdata)){
    if(param_realdata=="zeisel.imputed"){
      load('SymSim/match_params.zeisel_imputed.robj')
    }
    if(param_realdata=="pop4.romain"){
      load('SymSim/romain_pop4.zeisel.2.BP.robj')
      match_params <- cbind(result[[2]]$a[,1000,1],result[[2]]$b[,1000,1],result[[2]]$s[,1000,1])
    }
    match_params[,1]=log(base=10,match_params[,1])
    match_params[,2]=log(base=10,match_params[,2])
    match_params_den <- lapply(c(1:3),function(i){
      density(match_params[,i],n=2000)
    })
    if(joint==F){
      params <- Get_params(gene_effects,evf_res[[1]],match_params_den,bimod)      
    }else{params <- Get_params3(gene_effects,evf_res[[1]],match_params,bimod)}
  }else{
    params <- Get_params2(gene_effects,evf_res[[1]],bimod,ranges)
  }
  
  counts <- lapply(c(1:ngenes),function(i){
    count <- sapply(c(1:ncells_total),function(j){
      y <- rbeta(1,params[[1]][i,j],params[[2]][i,j])
      x <- rpois(1,y*params[[3]][i,j])
      return(x)
    })
  })
  cell_meta <- cbind( cellid=paste('cell',seq(1,ncells_total),sep='_'),evf_res[[2]],evf_res[[1]])
  counts <- do.call(rbind,counts)
<<<<<<< HEAD
  return(list(counts,gene_effects,evf_res,params))
=======
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
>>>>>>> update vignette files
=======
  if(SE==T){
    se <- SummarizedExperiment(assays=list(counts = as.matrix(counts),
                                           logcounts = log2(as.matrix(counts) + 1)),colData=cell_meta)
    rowData(se) <- gene_effects
    rowData(se)$gene_id <- paste('gene',seq(1,ngenes))
    return(se)
  }else{
    if (evf_type=="discrete")
      return(list(counts=counts,gene_effects=gene_effects,cell_meta=cell_meta,
                  params=params))
    else
      return(list(counts=counts,gene_effects=gene_effects,cell_meta=cell_meta,
                  params=params))
  }
}


<<<<<<< HEAD

#' Simulate observed count matrix given technical biases and the true counts
#' @param ncells_total number of cells
#' @param meta_cell the meta information related to cells, will be combined with technical cell level information and returned 
#' @param protocol a string, can be "ss2" or "umi"
#' @param alpha_mean the mean of rate of subsampling of transcripts during capture step, default at 10% efficiency
#' @param alpha_sd the std of rate of subsampling of transcripts
#' @param lenslope amount of length bias
#' @param nbins number of bins for gene length
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high, default is 0.8
#' @param nPCR the number of PCR cycles, default is 16
#' @param depth_mean mean of sequencing depth
#' @param depth_sd std of sequencing depth
#' @param SE input, should be a summerized experiment rather than a list of elements, default is False

True2ObservedCounts <- function(SE=NULL,true_counts,meta_cell,protocol,alpha_mean=0.1,alpha_sd=0.02,
                                lenslope=0.01,nbins=20,gene_len,amp_bias_limit=c(-0.2, 0.2),
                                rate_2PCR=0.8,nPCR=16,depth_mean, depth_sd){  
  if(!is.null(SE)){
    meta_cell <- colData(SE)
    true_counts <- assays(SE)$count
  }
  ngenes <- dim(true_counts)[1]; ncells <- dim(true_counts)[2]
  amp_bias <- cal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
  rate_2cap_vec <- rnorm(ncells, mean = alpha_mean, sd=alpha_sd)
  rate_2cap_vec[which(rate_2cap_vec < 0.0005)] <- 0.0005
  depth_vec <- rnorm(ncells, mean = depth_mean, sd=depth_sd)
  depth_vec[which(depth_vec < 200)] <- 200
  observed_counts <- matrix(0, ngenes, ncells)
  for (icell in 1:ncells){
    observed_counts[, icell] <- amplify_1cell(true_counts_1cell =  true_counts[, icell], protocol=protocol, 
                                              rate_2cap=rate_2cap_vec[icell], gene_len=gene_len, amp_bias = amp_bias, 
                                              rate_2PCR=rate_2PCR, nPCR=nPCR, N_molecules_SEQ = depth_vec[icell]) 
  }
  meta_cell2 <- data.frame(alpha=rate_2cap_vec,depth=depth_vec)
  meta_cell <- cbind(meta_cell, meta_cell2)
  if(is.null(SE)){return(list(observed_counts, meta_cell))}else{
    assays(SE)$observed_counts <- observed_counts
    colData(SE)<-meta_cell
    return(SE)
  }
>>>>>>> de68f0d77d0fb7c5402e989a978af800e9338ce1
}

#' Simulate observed count matrix given technical biases and the true counts in multiple batches
#' @param batch_effect the ratio of standard deviation to mean when sampling per batch parameters, default is 0.1
#' @param nbatch number of batches, default is 3
#' @param true_counts output of SimulateTrueCounts, with true count matrix and cell meta data
=======
# Batch_True2ObservedCounts <- function(){

<<<<<<< 176a52e550dca83672076caa19801f69fc1b3d4b
# }
# mean_alpha_mean <- rnorm()
# mean_lenslope <- rnorm()
#   ....
# lapply(c(1:nbatch),function(i){
#   True2ObservedCounts(alpha_mean=mean_alpha_mean[i],)
#   })
=======
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
>>>>>>> update vignette files

#' Simulate observed count matrix given technical biases and the true counts
#' @param ncells_total number of cells
#' @param meta_cell the meta information related to cells, will be combined with technical cell level information and returned 
#' @param nbatches number of batches (so far only 1)
>>>>>>> b1c727f6ac83412053a415eccd221e3a5179bbd7
#' @param protocol a string, can be "ss2" or "umi"
#' @param alpha_mean the mean of rate of subsampling of transcripts during capture step, default at 10% efficiency
#' @param alpha_sd the std of rate of subsampling of transcripts
#' @param lenslope amount of length bias
#' @param nbins number of bins for gene length
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high, default is 0.8
#' @param nPCR the number of PCR cycles, default is 16
#' @param depth_mean mean of sequencing depth
#' @param depth_sd std of sequencing depth
<<<<<<< HEAD
BatchTrue2ObservedCounts <- function(
  true_counts,protocol,batch_effect=0.1,nbatch=3,
  alpha_mean=0.1,alpha_sd=0.02,lenslope=0.01,
  nbins=20,amp_bias_limit=c(-0.2, 0.2),
  rate_2PCR=0.8,nPCR=18,depth_mean, depth_sd){
  ncells=length(true_counts[[1]][1,])
  ngenes = length(true_counts[[1]][,1])
  load("SymSim/gene_len_pool.RData")
  gene_len <- sample(gene_len_pool[which(gene_len_pool>100)], ngenes, replace = FALSE)
  batch_alpha_mean=rnorm(3,alpha_mean,alpha_mean*batch_effect)
  batch_lenslope=rnorm(3,lenslope,lenslope*batch_effect) 
  batch_rate_2PCR=rnorm(3,rate_2PCR,rate_2PCR*batch_effect) 
  batch_rate_2PCR[batch_rate_2PCR>1]=0.95
  batch_depth_mean=rnorm(3,depth_mean,depth_mean*batch_effect)
  if(sum(c(batch_alpha_mean,batch_rate_2PCR,batch_depth_mean)<0)){
    stop("batch effect is too large and created negative values in technical parameters")
  }
  ncells_batch <- rbinom(nbatch,p=1/nbatch,size=ncells)
  ncells_batch[3] <- ncells_batch[3]-(sum(ncells_batch)-ncells)
  batch_id <- do.call(c,lapply(c(1:nbatch),function(i){rep(i,ncells_batch[i])}))
  batch_id <- sample(batch_id,length(batch_id))
  counts <- lapply(c(1:nbatch),function(i){
    true_counts[[1]][,batch_id==i]
  })
  meta <- lapply(c(1:nbatch),function(i){
    true_counts[[3]][batch_id==i,]
  })
  observed_batches <- lapply(c(1:nbatch),function(i){
    observed_counts <- True2ObservedCounts(
      true_counts=counts[[i]],meta_cell=meta[[i]],
      protocol="umi",alpha_mean=batch_alpha_mean[i],lenslope=batch_lenslope[i],
      gene_len=gene_len,amp_bias_limit=c(-0.2, 0.2),
      rate_2PCR=batch_rate_2PCR[i],nPCR=18,
      depth_mean=batch_depth_mean[i], depth_sd=depth_sd)
    return(observed_counts)
  })
  batch_counts <- do.call(cbind,lapply(observed_batches,function(X){X[[1]]}))
  meta_cell <- do.call(rbind,lapply(observed_batches,function(X){X[[2]]}))
  meta_cell <- cbind(meta_cell,batch_id)
  observed_counts=batch_counts
  return(list(observed_counts,meta_cell))
}
=======
#' @param SE input, should be a summerized experiment rather than a list of elements, default is False
>>>>>>> b1c727f6ac83412053a415eccd221e3a5179bbd7

True2ObservedCounts <- function(SE=NULL,true_counts,meta_cell,nbatch=1,protocol,alpha_mean=0.1,alpha_sd=0.02,
                                lenslope=0.01,nbins=20,gene_len,amp_bias_limit=c(-0.2, 0.2),
                                rate_2PCR=0.8,nPCR=16,depth_mean, depth_sd){  
  if(!is.null(SE)){
    meta_cell <- colData(SE)
    true_counts <- assays(SE)$count
  }
  ngenes <- dim(true_counts)[1]; ncells <- dim(true_counts)[2]
  amp_bias <- cal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
  rate_2cap_vec <- rnorm(ncells, mean = alpha_mean, sd=alpha_sd)
  rate_2cap_vec[which(rate_2cap_vec < 0.0005)] <- 0.0005
  depth_vec <- rnorm(ncells, mean = depth_mean, sd=depth_sd)
  depth_vec[which(depth_vec < 200)] <- 200
  observed_counts <- matrix(0, ngenes, ncells)
  for (icell in 1:ncells){
    observed_counts[, icell] <- amplify_1cell(true_counts_1cell =  true_counts[, icell], protocol=protocol, 
                                              rate_2cap=rate_2cap_vec[icell], gene_len=gene_len, amp_bias = amp_bias, 
                                              rate_2PCR=rate_2PCR, nPCR=nPCR, N_molecules_SEQ = depth_vec[icell]) 
  }
  meta_cell2 <- data.frame(alpha=rate_2cap_vec,depth=depth_vec)
  meta_cell <- cbind(meta_cell, meta_cell2)
  if(is.null(SE)){return(list(observed_counts, meta_cell))}else{
    assays(SE)$observed_counts <- observed_counts
    colData(SE)<-meta_cell
    return(SE)
  }
}

#' Simulate technical biases 
#' @param lenslope amount of length bias
#' @param nbins number of bins for gene length
#' @param gene_len transcript length of each gene
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
cal_amp_bias <- function(lenslope, nbins, gene_len, amp_bias_limit){
  
  ngenes <- length(gene_len)
  len_bias_bin <- (-c(1:nbins))*lenslope
  len_bias_bin <- len_bias_bin-median(len_bias_bin)
  if (max(len_bias_bin) > amp_bias_limit[2]) {
    stop("The lenslope parameter is too large.")
  }
  max_rand_bias <- amp_bias_limit[2] - max(len_bias_bin)
  
  rand_bias <- rnorm(ngenes, mean=0, sd=max_rand_bias)
  rand_bias[rand_bias > max_rand_bias] <- max_rand_bias
  rand_bias[rand_bias < -max_rand_bias] <- -max_rand_bias
  #rand_bias <- runif(ngenes, -max_rand_bias,  max_rand_bias)
  
  binsize <- floor(ngenes/nbins)
  genes_in_bins <- vector("list", nbins)
  bin4genes <- numeric(ngenes)
  for (ibin in 1:(nbins-1)){
    genes_in_bins[[ibin]] <- order(gene_len)[((ibin-1)*binsize+1) : (ibin*binsize)]
    bin4genes[genes_in_bins[[ibin]]] <- ibin
  }
  genes_in_bins[[nbins]] <- order(gene_len)[((nbins-1)*binsize+1) : ngenes]
  bin4genes[genes_in_bins[[nbins]]] <- nbins
  
  len_bias <- numeric(ngenes); len_bias <- len_bias_bin[bin4genes]
  amp_bias <- rand_bias+len_bias
  return(amp_bias)
}

#' expand transcript counts to a vector of binaries of the same length of as the number of transcripts
#' @param true_counts_1cell number of transcript in one cell
expand2binary <- function(true_counts_1cell){
  expanded_vec <- rep(1, sum(true_counts_1cell))
  trans_idx <- sapply(which(true_counts_1cell>0), 
                      function(igene){return(rep(igene, true_counts_1cell[igene]))})
  trans_idx <- unlist(trans_idx)
  return(list(expanded_vec, trans_idx))
}
