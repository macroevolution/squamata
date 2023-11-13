basedir <- '~/Dropbox/Oz_Crown_Ages/dataArchive'

setwd(basedir)

setwd('./scripts/4.comparative_analyses/model-fits')

MAX <- 5


# Data: elongation, mass, vertebral 

source("1.load-data.R")


#--------------# 
# 1. elongation

D <- dataset_elongation(fulltree, basedir)
 
possible_shifts <- find_clades_pairs(D$phy, fulltree , clademat)
 
dmat <- getDescendantMatrix(D$phy)
 
res_el <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=MAX, ignore.cov=T, dmat=dmat)

summ_el <- summarizeStepwiseFits(res_el)	
write.table(summ_el, file = "fit_elong.csv", sep=",") 

#--------------# 
# 2. mass

D <- dataset_mass(fulltree, basedir)
 
possible_shifts <- find_clades_pairs(D$phy, fulltree , clademat)
 
dmat <- getDescendantMatrix(D$phy)
 
res_mass <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=MAX, ignore.cov=T, dmat=dmat)

summ_mass <- summarizeStepwiseFits(res_mass)
write.table(summ_mass, file = "fit_mass.csv", sep=",") 

# ------------- # 
# 3. vertebral

D <- dataset_vertebral(fulltree, basedir)

possible_shifts <- find_clades_pairs(D$phy, fulltree , clademat)
 
dmat <- getDescendantMatrix(D$phy)
 
res_vert <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=MAX, ignore.cov=T, dmat=dmat)

summ_vert <- summarizeStepwiseFits(res_vert)
write.table(summ_vert, file = "fit_vert.csv", sep=",") 


 