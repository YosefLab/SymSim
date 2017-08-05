required_packages <- c('roxygen2','devtools','FNN','reshape','ggplot2','gridExtra','ape','MASS','Rtsne','RColorBrewer','fAsianOptions','fOptions','fBasics','timeSeries','timeDate','ape' ,'reshape2','labeling')

install_packages <- rownames(installed.packages())
need_install <- setdiff(required_packages,install_packages)
if(length(need_install)>0){install.packages(need_install)}

lapply(required_packages, require, character.only = TRUE)

load('allsim.20170706.robj')
load('match_params.robj')


load_all('SCsimR')
ngenes=5000

simname <- '20170804.1batch.'

phyla<-rcoal(5)
jpeg(paste(simname,'tree.jpeg',sep=''))
plot(phyla)
dev.off()
cor_evf_mean<-vcv.phylo(phyla,cor=T)
varplot<-melt(cor_evf_mean)
varplot$V1<- factor( varplot$X1, levels = phyla$tip.label[order(cor_evf_mean[,1])])
varplot$V2<- factor( varplot$X2, levels = phyla$tip.label[order(cor_evf_mean[1,])])
p<- ggplot(data=varplot,aes(x=V1, y=V2, fill=value)) + geom_tile() + theme_bw()
ggsave(p,file = paste(simname,'varcovar.jpeg',sep='') )

gene_effects <- GeneEffects(ngenes=ngenes,nevf=10,randseed=0,sd=1,prob=0.3)


lapply(seq(0,1,0.1),function(bimod){
	lapply(seq(0.1,0.5,0.1),function(alpha){
		lapply(seq(0.1,1,0.1),function(evf_sd){
			lapply(seq(0.1,1,0.1),function(batch){
				result <- Npop1Batch(phyla=phyla,nevf=10,
					evf_sd=evf_sd, ncells=rep(500,5),
					randseed=0,gene_effects=gene_effects,bimod=bimod,alpha=alpha,
					alpha_sd=0.05,nbins=10,gcbias=0.1,lenbias=0.1,batch=batch, noise=0.1)
				save(result,file=paste(simname, '.bimod=',bimod,'.alpha=',alpha,'.evf_sd=',evf_sd,'.batch=',batch,'.robj',sep=''))
			})
		})
	})
})


