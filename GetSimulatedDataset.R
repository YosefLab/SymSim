library('devtools')
library('ape')
library('MASS')
load_all('SCsimR')
load('allsim.20170706.robj')
load('match_params.robj')
source('ContinuousTree.R')

ncells <- 1000
ngenes <- 20000
nbatch <- 3
nevf <-10
phyla <- Phyla5()

simname <- '20170905.allgenes'
for(beta in seq(0.5,1,0.1)){
	for (alpha in seq(0.01,0.15,0.01)){
		for (Sigma in seq(0.1,1,0.1)){
			for(batch in seq(0.1,1,0.1)){
				result=SimulateDataset(ncells,ngenes,nbatch,nevf,phyla=phyla,evf_type='Discrete',
					beta=beta,alpha=alpha,batch_sd=batch,Sigma=Sigma,gene_effect_prob=0.1)
				save(result,file=paste(simname, '.beta=',beta,'.alpha=',alpha,'.evf_sd=',Sigma,'.batch=',batch,'.discrete.robj',sep=''))
			}
		}	
	}
}


simname <- '20170905.allgenes'
for(beta in seq(0.5,1,0.1)){
	for (alpha in seq(0.01,0.15,0.01)){
		for (Sigma in seq(0.1,1,0.1)){
			for(batch in seq(0.1,1,0.1)){
				result=SimulateDataset(ncells,ngenes,nbatch,nevf,phyla=phyla,evf_type='Continuous',
					beta=beta,alpha=alpha,batch_sd=batch,Sigma=Sigma,gene_effect_prob=0.1)
				save(result,file=paste(simname, '.beta=',beta,'.alpha=',alpha,'.evf_sd=',Sigma,'.batch=',batch,'.continuous.robj',sep=''))
			}
		}	
	}
}

