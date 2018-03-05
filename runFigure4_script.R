sink("runFigure4.sh")
counter=0
for(seed in c(1:3)){
	for(Sigma in seq(0.05,0.2,0.05)){
		for(alpha in seq(0.05,0.2,0.05)){
			# for(alpha_sd in seq(0.01,0.1,0.03)){
				for(depth in seq(1e4,5e4,2e4)){
					# for(depth_sd in seq(5e3,2e4,5e3)){
						for(cv in c(1:3)){
						for(fold in seq(5,15,5)){
							counter=counter+1
							string <- paste('Rscript ~/Desktop/figure4_automated.R',seed, Sigma, alpha, alpha/cv, depth, depth/cv, alpha/fold, alpha_sd/fold, depth/fold, depth_sd/fold,counter)
							cat(string)
							cat("\n")
						}

					}
				}
			}
		}
	}

sink()
