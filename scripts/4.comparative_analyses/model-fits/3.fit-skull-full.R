

# 3 flavors of analysis:
#
#    (1) Find subsets of variables where we can fit the full model with covariances 
#       being estimated. This is preferred.
#   
#    (2) Repeat, but on principal component scores for species. This will at least
#       reduce some of the collinearity effect. Retain PC axes that account for 
#       90% of the variability in the data
#
#    (3) full analysis with covariances set to zero
#       This will absolutely inflate evidence for a rate shift on any
#        branch, by assuming nonindependence among variables that are highly
#        collinear
# 
 


source("1.load-data.R")

D <- dataset_skull(fulltree)


# Possible shifts in this tree:
#    these locations are constrained to 
#    those in the subfamily-level tree.
#    Thus, to reduce search space, no shifts *within* subfamilies

possible_shifts <- find_clades_pairs(D$phy, fulltree , clademat)
 

dmat <- getDescendantMatrix(D$phy)
 
 
#---------------------------# 
# Test:
# Random subsets

# Subsets approach 

test <- F
if (test)
{
raw_var <- ncol(D$y)
NVAR <- 10

ss <- sort(sample(ncol(D$y), NVAR))
D_sub <- D
D_sub$y <- D$y[,ss]

ff <- fit_shift(D_sub$y, D_sub$phy, ignore.cov=F)

res <- shift_search_step(D_sub$y, D_sub$phy, possible_shifts, KMAX=5, ignore.cov=F, dmat=dmat)

tmp <- summarizeStepwiseFits(res)
}

 
#

NSIMS <- 50

NVAR <- 10

KMAX <- 5

mm <- rep(NA, NSIMS)

simvec <- (2*KMAX + 1) * NSIMS

for (ii in 1:NSIMS){
	cat(ii, "\n")
	ss <- sort(sample(ncol(D$y), NVAR))
	D_sub <- D
	D_sub$y <- D$y[,ss]
 	res <- shift_search_step(D_sub$y, D_sub$phy, possible_shifts, KMAX=KMAX, ignore.cov=F, dmat=dmat)

	tmp <- summarizeStepwiseFits(res)
 	tmp <- data.frame(sim=rep(ii, nrow(tmp)), nvar=rep(NVAR, nrow(tmp)), tmp)
 	if (ii == 1){
 		dff <- tmp
 	}else{
 		dff <- rbind(dff, tmp)
 	}
 	write.table(dff, file = "fit_skull.csv", sep=",")
 	
}





