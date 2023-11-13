basedir <- '~/Dropbox/Oz_Crown_Ages/dataArchive'

setwd(basedir)

setwd('./scripts/4.comparative_analyses/model-fits')


#type <- "censored-only" 
type <- "all"
 
# ------------------------------- #
 
source("1.load-data.R")

D <- dataset_diet(fulltree, basedir)


# Possible shifts in this tree:
#    these locations are constrained to 
#    those in the subfamily-level tree.
#    Thus, to reduce search space, no shifts *within* subfamilies

possible_shifts <- find_clades_pairs(D$phy, fulltree, clademat)
 

dmat <- getDescendantMatrix(D$phy)
 
 
if (type != "censored-only")
{ 

res <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=5, ignore.cov=F, dmat=dmat)

summ_diet <- summarizeStepwiseFits(res)
 
write.table(summ_diet, file = "fit_diet.csv", sep=",") 

}


if (type == "censored-only")
{
	
KMAX <- 15
res <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=KMAX, ignore.cov=F, dmat=dmat, run.shift=F)
	

llvec <- rep(NA, KMAX+1)
aicvec <- llvec

llvec[1] <- res[[1]]$logL
aicvec[1] <- res[[1]]$aic

for (ii in 1:KMAX){
	llvec[ii+1] <- res[[3]][[ii]]$logL
	aicvec[ii+1] <- res[[3]][[ii]]$aic
}

}






