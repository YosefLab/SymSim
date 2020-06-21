
#############################################################
# Master Equation Related Functions
#############################################################

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
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params <- function(gene_effects,evf,match_param_den,bimod,scale_s){
  params <- lapply(1:3, function(iparam){evf[[iparam]] %*% t(gene_effects[[iparam]])})
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
  scaled_params[[3]] <- t(apply(scaled_params[[3]],2,function(x){x<-10^x}))*scale_s
  
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
}

#' @export
get_prob <- function(glength){
  if (glength >= 1000){prob <- 0.7} else{
    if (glength >= 100 & glength < 1000){prob <- 0.78}
    else if (glength < 100) {prob <- 0}
  }
  return(prob)
}

#' This function simulates the amplification, library prep, and the sequencing processes.
#' @param true_counts_1cell the true transcript counts for one cell (one vector)
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param rate_2cap the capture efficiency for this cell
#' @param gene_len gene lengths for the genes/transcripts, sampled from real human transcript length
#' @param amp_bias amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high
#' @param nPCR1 the number of PCR cycles
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param N_molecules_SEQ number of molecules sent for sequencing; sequencing depth
#' @return read counts (if protocol="nonUMI") or UMI counts (if protocol="UMI)
#' @export
amplify_1cell <- function(true_counts_1cell, protocol, rate_2cap, gene_len, amp_bias, 
                          rate_2PCR, nPCR1, nPCR2, LinearAmp, LinearAmp_coef, N_molecules_SEQ){
  ngenes <- length(gene_len)
  if (protocol=="nonUMI"){data(len2nfrag)} else 
    if(protocol=="UMI"){ } else
    {stop("protocol input should be nonUMI or UMI")}
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
  
  # pre-amplification:
  if (LinearAmp){
    PCRed_vec <- captured_vec*LinearAmp_coef
  } else {
    temp <- runif(length(captured_vec)) < amp_rate
    temp <- temp*2+captured_vec-temp
    for (iPCR in 2:nPCR1){
      eff <- runif(length(temp))*amp_rate
      v1 <- temp*(1-eff)
      round_down <- (v1-floor(v1)) < runif(length(v1))
      v1[round_down] <- floor(v1[round_down]); v1[!round_down] <- ceiling(v1[!round_down])
      temp <- v1 + 2*(temp-v1)
    }
    PCRed_vec <- temp
  }
  
  if (protocol=="nonUMI"){ # add fragmentation step here
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
    # another 8 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2){
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
    }
    for (iPCR in 3:nPCR2){
      frag_vec <- frag_vec + round(frag_vec*rate_2PCR)
    }
    SEQ_efficiency=N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1) {read_count <- frag_vec} else{
      read_count <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)}) }
    return(read_count)
  } else if (protocol=="UMI"){
    
    prob_vec <- sapply(gene_len[trans_idx[1:(length(trans_idx)-1)]], get_prob)
    # fragmentation: 
    frag_vec <- sapply(1:(length(PCRed_vec)-1), function(igene)
    {return(rbinom(n=1, size = PCRed_vec[igene], prob = prob_vec[igene] ))})
    
    # another 10 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2){
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
    }
    
    frag_vec <- round(frag_vec * (1+rate_2PCR)^(nPCR2-1))
    
    SEQ_efficiency <- N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1){sequenced_vec <- frag_vec} else {
      sequenced_vec <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)})}
    
    temp_vec <- c(sequenced_vec,1)
    for (i in seq(2,1,-1)){
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec; 
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec)-1)]
    
    UMI_counts=numeric(ngenes); 
    GI=c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell>0)){
      x=recovered_vec[(GI[i]+1):GI[i+1]];
      UMI_counts[i]=sum(x>0); 
    }
    
    return(list(UMI_counts, sequenced_vec, sum(frag_vec>0)))
  }
}

#' sample from truncated normal distribution
#' @param a the minimum value allowed 
#' @param b the maximum value allowed
rnorm_trunc <- function(n, mean, sd, a, b){
  vec1 <- rnorm(n, mean = mean, sd=sd)
  beyond_idx <- which(vec1 < a | vec1 > b)
  if (length(beyond_idx) > 0) { # for each value < rate_2cap_lb
    substi_vec <- sapply(1:length(beyond_idx), function(i){
      while (TRUE){
        temp <- rnorm(1, mean = mean, sd=sd)
        if (temp > a & temp < b) {break}}
      return(temp)} )
    vec1[beyond_idx] <- substi_vec
  }
  return(vec1)
}

#' Compute value of impulse function given parameters.
#' Enforces lower bound on value of function to avoid numerical
#' errors during model fitting.
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


#' Creating an example tree with 5 tips
#' @param plotting True for plotting the tree on console, False for no plot 
#' @return a tree object
#' @export
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

#' Creating an example tree with 3 tips
#' @param plotting True for plotting the tree on console, False for no plot 
#' @return a tree object
#' @export
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
#' @param impulse if the impluse model should be used instead of Brownian motion
#' @param evf_center the mean of Gaussain function where the non-Diff EVFs are sampled from
#' @param vary which parameters are affected by Diff-EVFs. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s"
#' @param tip The leaf that the path with impulse lead to
#' @param Sigma The standard deviation of the brownian motion of EVFs changing along the tree 
#' @param plotting Whether to plot the trajectory or not
#' @param plotname The string to be used in the output file name
#' @param seed the random seed 
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the branch each cell comes from (pop) and its depth in the tree (depth)
ContinuousEVF <- function(phyla,ncells,n_nd_evf,n_de_evf,impulse=F,evf_center=1,vary='s',
                          Sigma,plotting=T,plotname='cont_evf.pdf',seed){
  set.seed(seed)
  edges <- cbind(phyla$edge,phyla$edge.length)
  edges <- cbind(c(1:length(edges[,1])),edges)
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
  neutral <- SampleSubtree(root,0,evf_center,edges,ncells)
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
          SampleSubtree(root,0,evf_center,edges,ncells,neutral=neutral)    
        })
        # pdf(file = plotname,width=15,height=5)
        # if(plotting==T){PlotRoot2Leave(cbind(neutral,tips,edges,root,internal)}
        # dev.off()
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
#' @param n_nd_evf number of non-Diff EVFs
#' @param n_de_evf number of Diff EVFs
#' @param vary which parameters are affected by Diff-EVFs. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s"
#' @param evf_center the value used to generated evf means. Suggested value is 1
#' @param seed the random seed
#' @return a list of two object, one is the evf, and the other is a dataframe indicating the population each cell comes from (pop)
DiscreteEVF <- function(phyla, ncells_total, min_popsize, i_minpop, Sigma, n_nd_evf, n_de_evf, 
                        vary, evf_center, seed){
  set.seed(seed)
  npop <- length(phyla$tip.label)
  # set the number of cells in each population: first give each population min_popsize cells
  # then randomly distribute the rest of cells to all populations except the smallest one
  ncells_pop <- rep(min_popsize, npop)
  if (ncells_total < min_popsize*npop) {
    stop("The size of the smallest population is too big for the total number of cells")}
  larger_pops <- setdiff(1:npop, i_minpop)
  ncells_pop[larger_pops] <- floor((ncells_total-min_popsize)/length(larger_pops))
  leftover <- ncells_total-sum(ncells_pop)
  if (leftover > 0){
    temp <- sample(larger_pops, leftover, replace = F); ncells_pop[temp] <- ncells_pop[temp] + 1
  }
  
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
}

#' Generate both evf and gene effect and simulate true transcript counts
#' @param ncells_total number of cells
#' @param min_popsize the number of cells in the smallest population
#' @param i_minpop specifies which population has the smallest size
#' @param ngenes number of genes
#' @param evf_center the value which evf mean is generated from
#' @param nevf number of evfs
#' @param evf_type string that is one of the following: 'one.population','discrete','continuous'
#' @param n_de_evf number of differential evfs between populations
#' @param vary which kinetic parameters should the differential evfs affect. Default is 's'. Can be "kon", "koff", "s", "all", "except_kon", "except_koff", "except_s". Suggestions are "all" or "s".
#' @param impulse use the impulse function when generating continuous population or not. Default is F. 
#' @param Sigma parameter of the std of evf values within the same population
#' @param phyla the cell developmental tree if chosing 'discrete' or 'continuous' evf type. Can either be generated randomly (using pbtree(nclusters) function from phytools package) or read from newick format file using the ape package
#' @param param_realdata pick from zeisel.imputed or NULL; zeisel.imputed means using the distribution of kinetic parameters learned from the Zeisel 2015 dataset. This option is recommended.
#' @param gene_effect_prob the probability that the effect size is not 0
#' @param geffect_mean the mean of the normal distribution where the non-zero gene effect sizes are sampled from 
#' @param gene_effect_sd the standard deviation of the normal distribution where the non-zero gene effect sizes are sampled from 
#' @param bimod the amount of increased bimodality in the transcript distribution, 0 being not changed from the results calculated using evf and gene effects, and 1 being all genes are bimodal
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @param prop_hge the proportion of very highly expressed genes
#' @param mean_hge the parameter to amplify the gene-expression levels of the very highly expressed genes
#' @param randseed random seed
#' @return a list of 4 elements, the first element is true counts, second is the gene level meta information, the third is cell level meta information, including a matrix of evf and a vector of cell identity, and the fourth is the parameters kon, koff and s used to simulation the true counts
#' @import phytools
#' @export
SimulateTrueCounts <- function(ncells_total,min_popsize,i_minpop=1,ngenes, 
                               evf_center=1,evf_type="one.population",nevf=10,
                               phyla, randseed, n_de_evf=0,vary='s',Sigma=0.4,
                               geffect_mean=0,gene_effects_sd=1,gene_effect_prob=0.3,
                               bimod=0,param_realdata="zeisel.imputed",scale_s=1,impulse=F,
                               prop_hge=0.015, mean_hge=5){
  set.seed(randseed)
  n_nd_evf=nevf-n_de_evf
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
    evf_res <- ContinuousEVF(phyla,ncells_total,n_nd_evf=nevf-n_de_evf,n_de_evf=n_de_evf,
                             evf_center=evf_center,vary=vary,impulse=impulse,
                             Sigma,plotting=T,seed=seed[1])    
  }
  gene_effects <- GeneEffects(ngenes=ngenes,nevf=nevf,randseed=seed[2],prob=gene_effect_prob,
                              geffect_mean=geffect_mean,geffect_sd=gene_effects_sd)
  if(!is.null(param_realdata)){
    if(param_realdata=="zeisel.imputed"){
      data(param_realdata.zeisel.imputed)
    } else {stop("wrong input for parameter param_realdata")}
    
    match_params[,1]=log(base=10,match_params[,1])
    match_params[,2]=log(base=10,match_params[,2])
    match_params[,3]=log(base=10,match_params[,3])
    match_params_den <- lapply(c(1:3),function(i){
      density(match_params[,i],n=2000)
    })
    params <- Get_params(gene_effects,evf_res[[1]],match_params_den,bimod,scale_s=scale_s)
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
  
  if (prop_hge>0){
    chosen_hge <- sample(ngenes, ceiling(ngenes*prop_hge), replace = F)
    multi_factors <- numeric(length(chosen_hge))
    rank_sum <- rank(rowSums(params[[3]][chosen_hge,]))
    multi_factors <- sapply(1:length(chosen_hge), function(igene){
      tosubstract <- -rank_sum[igene]*1/length(chosen_hge)+1
      if (runif(1, 0, 1) < 1){
        multi_factor <- mean_hge - tosubstract
      } else {multi_factor <- mean_hge}
      return(multi_factor)
    })
    new_s <- matrix(0, length(chosen_hge), ncells_total)
    for (i in 1:length(chosen_hge)){
      new_s[i,] <- params[[3]][chosen_hge[i],] * (2^multi_factors[i])
    }
    params[[3]][chosen_hge,] <- new_s
    counts[chosen_hge] <- lapply(1:length(chosen_hge),function(i){
      s_vec <- new_s[i,]
      count <- sapply(c(1:ncells_total),function(j){
        x <- rpois(1,s_vec[j])
        return(x)
      })
      return(count)
    })
    chosen_hge <- cbind(chosen_hge, multi_factors)
  } else {chosen_hge <- NULL}
  
  cell_meta <- cbind( cellid=paste('cell',seq(1,ncells_total),sep='_'),evf_res[[2]],evf_res[[1]])
  counts <- do.call(rbind,counts)
  
  return(list(counts=counts,gene_effects=gene_effects,cell_meta=cell_meta,kinetic_params=params))
}


#' Simulate observed count matrix given technical biases and the true counts
#' @param true_counts gene cell matrix
#' @param meta_cell the meta information related to cells, will be combined with technical cell level information and returned 
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param alpha_mean the mean of rate of subsampling of transcripts during capture step, default at 10 percent efficiency
#' @param alpha_sd the std of rate of subsampling of transcripts
#' @param lenslope amount of length bias
#' @param nbins number of bins for gene length
#' @param gene_len a vector with lengths of all genes
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high, default is 0.8
#' @param nPCR1 the number of PCR cycles in "pre-amplification" step, default is 16
#' @param nPCR2 the number of PCR cycles used after fragmentation. 
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param depth_mean mean of sequencing depth
#' @param depth_sd std of sequencing depth
#' @import SummarizedExperiment
#' @export
True2ObservedCounts <- function(true_counts,meta_cell,protocol,alpha_mean=0.1,alpha_sd=0.002,
                                gene_len,depth_mean, depth_sd, lenslope=0.02,nbins=20,
                                amp_bias_limit=c(-0.2, 0.2),
                                rate_2PCR=0.8,nPCR1=16, nPCR2=10, LinearAmp=F, LinearAmp_coef=2000){  
  # if(!is.null(SE)){
  #   meta_cell <- colData(SE)
  #   true_counts <- assays(SE)$count
  # }
  ngenes <- dim(true_counts)[1]; ncells <- dim(true_counts)[2]
  amp_bias <- cal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
  rate_2cap_lb <- 0.0005; depth_lb <- 200 # lower bound for capture efficiency and sequencing depth  
  rate_2cap_vec <- rnorm_trunc(n=ncells, mean = alpha_mean, sd=alpha_sd, a=rate_2cap_lb, b=1)
  depth_vec <- rnorm_trunc(n=ncells, mean = depth_mean, sd=depth_sd,a=depth_lb,b=Inf)

  observed_counts <- lapply(c(1:ncells),function(icell){
    amplify_1cell(true_counts_1cell =  true_counts[, icell], protocol=protocol, 
                  rate_2cap=rate_2cap_vec[icell], gene_len=gene_len, amp_bias = amp_bias, 
                  rate_2PCR=rate_2PCR, nPCR1=nPCR1, nPCR2=nPCR2, LinearAmp = LinearAmp, 
                  LinearAmp_coef = LinearAmp_coef, N_molecules_SEQ = depth_vec[icell])     
  })

  meta_cell2 <- data.frame(alpha=rate_2cap_vec,depth=depth_vec,stringsAsFactors = F)
  meta_cell <- cbind(meta_cell, meta_cell2)
  
  if (protocol=="UMI"){
    UMI_counts <- do.call(cbind, lapply(observed_counts, "[[", 1))
    nreads_perUMI <- lapply(observed_counts, "[[", 2)
    nUMI2seq <- sapply(observed_counts, "[[", 3)
    observed_counts <- UMI_counts
  } else
    observed_counts <- do.call(cbind,observed_counts)
  
  if (protocol=="UMI"){return(list(counts=observed_counts, cell_meta=meta_cell, nreads_perUMI=nreads_perUMI, 
                                   nUMI2seq=nUMI2seq))
  } else
    return(list(counts=observed_counts, cell_meta=meta_cell))
}

#' Divide the observed counts into multiple batches by adding batch effect to each batch
#' @param observed_counts_res the output from True2ObservedCounts
#' @param nbatch number of batches
#' @param batch_effect_size amount of batch effects. Larger values result in bigger differences between batches. Default is 1.
#' @export
DivideBatches <- function(observed_counts_res, nbatch, batch_effect_size=1){
  ## add batch effects to observed counts
  # use different mean and same sd to generate the multiplicative factor for different gene in different batch
  observed_counts <- observed_counts_res[["counts"]]
  meta_cell <- observed_counts_res[["cell_meta"]]
  ncells <- dim(observed_counts)[2]; ngenes <- dim(observed_counts)[1]
  batchIDs <- sample(1:nbatch, ncells, replace = TRUE)
  meta_cell2 <- data.frame(batch=batchIDs, stringsAsFactors = F)
  meta_cell <- cbind(meta_cell, meta_cell2)
  
  mean_matrix <- matrix(0, ngenes, nbatch)
  gene_mean <- rnorm(ngenes, 0, 1)
  temp <- lapply(1:ngenes, function(igene) {
    return(runif(nbatch, min = gene_mean[igene]-batch_effect_size, max = gene_mean[igene]+batch_effect_size))
  })
  mean_matrix <- do.call(rbind, temp)
  
  batch_factor <- matrix(0, ngenes, ncells)
  for (igene in 1:ngenes){
    for (icell in 1:ncells){
      batch_factor[igene, icell] <- rnorm(n=1, mean=mean_matrix[igene, batchIDs[icell]], sd=0.01)
    }
  }
  observed_counts <- round(2^(log2(observed_counts)+batch_factor))
  return(list(counts=observed_counts, cell_meta=meta_cell))
}

#' Simulate technical biases 
#' @param lenslope amount of length bias. This value sould be less than 2*amp_bias_limit[2]/(nbins-1)
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
