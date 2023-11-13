basedir <- '~/Dropbox/Oz_Crown_Ages/dataArchive'

setwd(basedir)

setwd('./scripts/4.comparative_analyses/model-fits')

# type <- "censored-only" 
type <- "all"
 
# Model fitting for non-phylo diet  

# ------------------------------- #
 
source("1.load-data.R")

# ----------------------
# FULL DIET DATASET

D <- dataset_diet_noTree(fulltree, basedir)

# Possible shifts in this tree:
#    these locations are constrained to 
#    those in the subfamily-level tree.
#    Thus, to reduce search space, no shifts *within* subfamilies

possible_shifts <- find_clades_pairs(D$phy, fulltree, clademat)
 
dmat <- getDescendantMatrix(D$phy)
 
if (type != "censored-only")
{ 

res <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=3, ignore.cov=F, dmat=dmat)

summ_diet <- summarizeStepwiseFits(res)
 
write.table(summ_diet, file = "fit_diet_noTree.csv", sep=",") 

}


# if (type == "censored-only")
# {
	
# KMAX <- 5
# res <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=KMAX, ignore.cov=F, dmat=dmat, run.shift=F)
	

# llvec <- rep(NA, KMAX+1)
# aicvec <- llvec

# llvec[1] <- res[[1]]$logL
# aicvec[1] <- res[[1]]$aic

# for (ii in 1:KMAX){
	# llvec[ii+1] <- res[[3]][[ii]]$logL
	# aicvec[ii+1] <- res[[3]][[ii]]$aic
# }

# }



# -------------------------
# PCA OF FULL DIET DATASET

D <- dataset_dietPC_noTree(fulltree, basedir)

# Possible shifts in this tree:
#    these locations are constrained to 
#    those in the subfamily-level tree.
#    Thus, to reduce search space, no shifts *within* subfamilies

possible_shifts <- find_clades_pairs(D$phy, fulltree, clademat)
 

dmat <- getDescendantMatrix(D$phy)
 
 
if (type != "censored-only")
{ 

res <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=3, ignore.cov=F, dmat=dmat)

summ_diet <- summarizeStepwiseFits(res)
 
write.table(summ_diet, file = "fit_diet_PC_noTree.csv", sep=",") 

}





# ----------------------
# RANDOM SUBSETS


NSIMS <- 10

NVAR <- 15

KMAX <- 5

D <- dataset_diet_noTree(fulltree, basedir)
dmat <- getDescendantMatrix(D$phy)
possible_shifts <- find_clades_pairs(D$phy, fulltree, clademat)


for (ii in 1:NSIMS){
    cat(ii, "\n")
    ss <- sort(sample(ncol(D$y), NVAR))
    D_sub <- D
    D_sub$y <- D$y[,ss]
    res <- shift_search_step(D_sub$y, D_sub$phy, possible_shifts, KMAX=KMAX, ignore.cov=F, dmat=dmat, run.shift=T)

    tmp <- summarizeStepwiseFits(res)
    tmp <- data.frame(sim=rep(ii, nrow(tmp)), nvar=rep(NVAR, nrow(tmp)), tmp)
    if (ii == 1){
        dff <- tmp
    }else{
        dff <- rbind(dff, tmp)
    }
}

write.table(dff, file = "fit_diet_noTree_randomSubset.csv", sep=",")


