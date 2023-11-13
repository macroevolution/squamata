basedir <- '~/Dropbox/Oz_Crown_Ages/dataArchive'

setwd(basedir)

setwd('./scripts/4.comparative_analyses/model-fits')


source("1.load-data.R")

D <- dataset_climate(fulltree, basedir)


# Possible shifts in this tree:
#    these locations are constrained to 
#    those in the subfamily-level tree.
#    Thus, to reduce search space, no shifts *within* subfamilies

possible_shifts <- find_clades_pairs(D$phy, fulltree, clademat)
 

dmat <- getDescendantMatrix(D$phy)
 
res <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=10, ignore.cov=F, dmat=dmat)

summ_climate <- summarizeStepwiseFits(res)
 
write.table(summ_climate, file = "fit_climate.csv", sep=",") 





