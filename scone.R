library(scone)
library(RColorBrewer)
cc <- c(brewer.pal(9, "Set1"))

# load('alpha_change.robj')
# dropped<-do.call(cbind,lapply(alpha_change,function(X){X[[3]][,c(1:100)]}))
# true <-  do.call(cbind,lapply(alpha_change,function(X){X[[2]][,c(1:100)]}))
# normed <- Normalize(true)

sim1 <- temp
genenames = paste('gene',as.character(c(1:length(sim1[[2]][,1]))),sep='_')
samplename = paste('cell',as.character(c(1:length(sim1[[2]][,1]))),sep='_')
batch = factor(rep(1,length(samplename)))

for(i in c(2:4)){
  colnames(sim1[[i]]) = samplename
  rownames(sim1[[i]]) = genenames
}

## ----fnr_fit-------------------------------------------------------------
# Extract Housekeeping Genes

hk <- genenames[rank(rowMeans(sim1[[2]]))>900]
# Mean log10(x+1) expression
counts <- sim1[[4]]
mu_obs = rowMeans(log10(counts[hk,]+1))

# Assumed False Negatives
drop_outs = (counts[hk,] == 0)

# Logistic Regression Model of Failure
ref.glms = list()
for (si in 1:dim(drop_outs)[2]){
  fit = glm(cbind(drop_outs[,si],1 - drop_outs[,si]) ~ mu_obs,
            family=binomial(logit))
  ref.glms[[si]] = fit$coefficients
}


## ----fnr_vis,fig.width=8,fig.height=4,out.width="800px",out.height="400px"----

par(mfrow=c(1,2))

# Plot Failure Curves and Calculate AUC
plot(NULL, main = "False Negative Rate Curves",
     ylim = c(0,1),xlim = c(0,6), 
     ylab = "Failure Probability", xlab = "Mean log10 Expression")
x = (0:60)/10
AUC = NULL
for(si in 1:20){
  y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
  AUC[si] = sum(y)/10
  lines(x, 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1),
        type = 'l', lwd = 2, col = cc[batch][si])
}

# Barplot of FNR AUC
o = order(AUC)[order(batch[order(AUC)])]

barplot(AUC[o], col=cc[batch][o], border=cc[batch][o], main="FNR AUC")
legend("topright", legend=levels(batch), fill=cc, cex=0.4)


## ----metric_sample_filter, fig.height= 10,out.height="1000px"------------

# Initial Gene Filtering: 
# Select "common" transcripts based on proportional criteria.
num_reads = quantile(counts[counts > 0])[4]
num_cells = 0.25*ncol(counts)
is_common = rowSums(counts >= num_reads ) >= num_cells

# Metric-based Filtering
mfilt = metric_sample_filter(counts,
                             gene_filter = is_common,
                             pos_controls = rownames(counts) %in% hk,
                             zcut = 3, mixture = FALSE,
                             plot = TRUE)

# Simplify to a single logical
mfilt = !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)

goodDat = counts[,mfilt]

# Final Gene Filtering: Highly expressed in at least 5 cells
num_reads = quantile(counts[counts > 0])[4]
num_cells = 5
is_quality = rowSums(counts >= num_reads ) >= num_cells

# Biological Origin - Variation to be preserved (Optional)
bio = factor(rep(1,length(samplename)))

# Creating a SconeExperiment Object
my_scone <- SconeExperiment(expr,
                qc=ppq, bio = bio,
                negcon_ruv = rownames(expr) %in% negcon,
                poscon = rownames(expr) %in% poscon
)
my_scone <- SconeExperiment(counts,
                # negcon_ruv = rownames(expr) %in% negcon,
                # poscon = rownames(expr) %in% poscon,qc=ppq, 
                bio = bio
)


## ----scone_in2-----------------------------------------------------------

## ----- User-defined function: Dividing by number of detected genes -----

EFF_FN = function (ei)
{
  sums = colSums(ei > 0)
  eo = t(t(ei)*sums/mean(sums))
  return(eo)
}

## ----- Scaling Argument -----

scaling=list(none=identity, # Identity - do nothing
             eff = EFF_FN, # User-defined function
             sum = SUM_FN, # SCONE library wrappers...
             tmm = TMM_FN, 
             uq = UQ_FN,
             fq = FQT_FN,
             deseq = DESEQ_FN)


## ----scone_in3, eval=FALSE-----------------------------------------------
#  
 # Simple FNR model estimation with SCONE::estimate_ziber
 fnr_out = estimate_ziber(x = counts, bulk_model = TRUE,
                          pos_controls = rownames(counts) %in% hk,
                          maxiter = 100)
 
 ## ----- Imputation List Argument -----
 imputation=list(none=impute_null, # No imputation
                 expect=impute_expectation) # Replace zeroes
 
 ## ----- Imputation Function Arguments -----
 # accessible by functions in imputation list argument
 impute_args = list(p_nodrop = fnr_out$p_nodrop, mu = exp(fnr_out$Alpha[1,]))
 
## ----scone_run-----------------------------------------------------------
BiocParallel::register(
  BiocParallel::SerialParam()
) # Register BiocParallel Serial Execution

my_scone <- scone(my_scone,
                 imputation = imputation, impute_args = impute_args,
                  scaling=scaling,
                  run=TRUE,
                  k_qc=0, k_ruv = 0,
                  eval_kclust = 2:6,stratified_pam = TRUE,
                  return_norm = "in_memory",
                 adjust_bio="no",
                  zero = "preadjust")




## ----scone_params_view---------------------------------------------------

apply(get_params(my_scone),2,unique)


## ----scone_view1---------------------------------------------------------

# View Metric Scores
head(get_scores(my_scone))

# View Mean Score Rank
head(get_score_ranks(my_scone))

# Extract normalized data from top method
out_norm <- lapply(c(1:length(rownames(get_params(my_scone)))),function(i){
  get_normalized(my_scone,method = rownames(get_params(my_scone))[i])
})


## ----biplot_color--------------------------------------------------------

pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                center = TRUE,scale = FALSE)
bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)


## ----biplot_color4-------------------------------------------------------

bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)

points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)

points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
       pch = 1, col = "blue", cex = 1)
points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
       pch = 1, col = "blue", cex = 1.5)

arrows(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][1],
       bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][2],
       bp_obj[1,][1],
       bp_obj[1,][2],
       lty = 2, lwd = 2)


## ----sconeReport, eval=FALSE---------------------------------------------
#  
#  # Methods to consider
 scone_methods = c(rownames(get_params(my_scone))[1:12],
                   "none,none,no_uv,no_bio,no_batch")
 
 # Shiny app
 sconeReport(my_scone,methods = scone_methods,
             qc = ppq,
             bio = bio,
             negcon = negcon, poscon = poscon)
 

## ----session-------------------------------------------------------------
sessionInfo()
