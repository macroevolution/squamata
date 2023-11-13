rm(list = ls())
# --------------------------# 
# Should be self-contained example file 
#   illustrating model-fitting on the skull dataset
#
#   TODO: Needs cleanup
#
#

# Worked example with skull dataset:

source("1.load-data.R")

# construct target dataset:
D <- dataset_skull(fulltree)
 
# Possible shifts in this tree:
#    these locations are constrained to 
#    those in the subfamily-level tree.
#    Thus, to reduce search space, no shifts *within* subfamilies

possible_shifts <- find_clades_pairs(D$phy, fulltree , clademat)
 

dmat <- getDescendantMatrix(D$phy)
 
 
res <- shift_search_step(D$y[,c(1,3,5, 7 , 9)], D$phy, possible_shifts, KMAX=5, ignore.cov=F, dmat=dmat)
 
tmp <- summarizeStepwiseFits(res)

#----------------------------# 

# example fitting 0 shift + 1 + 2 covariance regimes 
#   to skull data, with censored model 
#    (eg, shift in cov matrix + shift in trait means)
library(ape)
snakes <- getMRCA(D$phy, tip = c("Sonora_semiannulata", "Typhlops_jamaicensis"))
colub <- getMRCA(D$phy, tip = c("Sonora_semiannulata", "Xenodermus_javanicus"))

traits <- sample(1:40, size=10)
fit0 <- fit_shift(D$y[, traits], D$phy, ignore.cov=F)
fit1 <- fit_censored(D$y[, traits], D$phy, shift=snakes, ignore.cov=F)
fit2 <- fit_censored(D$y[, traits], D$phy, shift=c(snakes, colub), ignore.cov=F)



#-------------------- #

D_pc <- D
pp <- princomp(D$y)
D_pc$y <- pp$scores[,1:2] #90% of variance

res_pc <- shift_search_step(D_pc$y, D_pc$phy, possible_shifts, KMAX=2, ignore.cov=T, dmat=dmat) 
tmp_pc <- summarizeStepwiseFits(res_pc) 
cbind(tmp_pc[,1:5], tmp[,1:5])
# -------------------- 

# 3 flavors of analysis:
# 
#   (1) full analysis with covariances set to zero
#       This will absolutely inflate evidence for a rate shift on any
#        branch, by assuming nonindependence among variables that are highly
#        collinear
# 
#   (2) Repeat, but on principal component scores for species. This will at least
#       reduce some of the collinearity effect. Retain PC axes that account for 
#       90% of the variability in the data
#
#   (3) Find subsets of variables where we can fit the full model with covariances 
#       being estimated. This is preferred.




#-----------------------------------------# 
# Subsets approach
D_sub <- D
raw_var <- ncol(D$y)
NVAR <- 8


ss <- sample(1:raw_var, NVAR)
D_sub$y <- D$y[,ss]
res_sub <- shift_search_step(D_sub$y, D_sub$phy, possible_shifts, KMAX=4, ignore.cov=F, dmat=dmat) 
summarizeStepwiseFits(res_sub)




pdf(file = "tree.pdf", height=30, width=5)
plot.phylo(D$phy, cex=0.1)
nodelabels(node = 1321, pch=19, col="red", cex=3)
dev.off()
#
D <- dataset_diet(fulltree)

possible_shifts <- find_clades_pairs(D$phy, fulltree , clademat)
dmat <- getDescendantMatrix(D$phy)
res <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=2, ignore.cov=T, dmat=dmat)






fit0 <- fit_shift(D$y, D$phy)
shifts <- find_clades(D$phy, fulltree, "Serpentes")
fit1_C <- fit_censored(D$y, D$phy, shift=shifts, ignore.cov=TRUE)
fit1_NC <- fit_shift(D$y, D$phy, shift=shifts, ignore.cov=TRUE)

best_C <- fit1_C

source("shift-search.R")
res <- shift_search_step(D$y, D$phy, possible_shifts, KMAX=10)

ss <- getShiftNodes(res[[2]][[1]])
 
D$y <- D$y[ , 1:35]
fit1_NC <- fit_shift(D$y, D$phy, shift = ss, ignore.cov = F, auto.drop=T) 
fit0 <- fit_shift(D$y, D$phy, ignore.cov=F, auto.drop = T)

fit1_C <- fit_censored(D$y, D$phy, shift=ss, ignore.cov=F, auto.drop=T)


length(res[[1]])

aicvec <- numeric(11)
aicvec[1] <- res[[1]]$aic
for (ii in 1:10){
	aicvec[ii+1] <- res[[2]][[ii]]$aic
}

# Try a worked example to contrast impact of ignore.cov = T vs F

fit0 <- fit_shift(D$y, D$phy)
shifts <- find_clades(D$phy, fulltree, "Serpentes")

fit1_no_cov <- fit_shift(D$y, D$phy, shift=shifts, ignore.cov=TRUE)
fit1_cov    <- fit_shift(D$y, D$phy, shift = shifts, ignore.cov=FALSE)

fit1_c_no_cov <- fit_censored(D$y, D$phy, shifts, ignore.cov=TRUE)
fit1_c_cov    <- fit_censored(D$y, D$phy, shifts, ignore.cov = FALSE)


aicvec_c <- rep(NA, length(possible_shifts))
aicvec   <- rep(NA, length(possible_shifts))




fit_censored(D$y, D$phy, shift = possible_shifts[1])



for (ii in 1:length(possible_shifts)){
	cat(ii, "\n")
 	fit_censored(D$y, D$phy, shift = possible_shifts[1])
}





 
pshifts2 <- find_clades_pairs(D$phy, fulltree , clademat)

fit0 <- fit_shift(D$y, D$phy) #zero shift likelihood
llvec <- numeric(length(pshifts2))
llvec_c <- llvec

for (ii in 1:length(pshifts2)){
	cat(ii, "\n")
	#tmpfit <- fit_shift(D$y, D$phy, possible_shifts[ii])
	tmpfit <- fit_censored(D$y, D$phy, pshifts2[ii])
	llvec_c[ii] <- tmpfit$logL
	
}

 
 
llvec_c <- llvec_c - min(llvec_c) 
 
keep <- pshifts2[llvec_c > 100]
 

# maybe throw away all observations within 500 of the MIN?
 
cset <- combn(keep, 5)



#---------------------------# 
# Example with skull
shifts <- find_clades(D$phy, fulltree, c("Serpentes", "Macrostomata"))
fit0 <- fit_shift(D$y, D$phy)   # omit shifts arg for constant rate
fit1 <- fit_shift(D$y, D$phy, shifts)
fit2 <- fit_censored(D$y, D$phy, shifts)

plot(fit1$phy.scaled, show.tip.label=FALSE)
plot(fit2$phy.scaled, cex=0.5)

fit5 <- fit_censored(D$y, D$phy, cset[,1])

for (ii in 1:100){
	cat(ii, "\n")
	tmp <- fit_censored(D$y, D$phy, cset[,ii])
}


cset[,18]
shifts <- cset[,18]

fitx <- fit_censored(D$y, D$phy, shifts[3:5])

t1 <- ape::extract.clade(D$phy, node = shifts[3])
t2 <- ape::extract.clade(D$phy, node = shifts[4])
t3 <- ape::extract.clade(D$phy, node = shifts[5])

# 2 procedures: 
#   (1) Stepwise
#   (2) multi-model, not stepwise

plot(D$phy, show.tip.label=F)
ape::nodelabels(node = shifts[3:5], pch=19, col="red", cex=3)

dd <- getDescendantMatrix(fulltree)



fit0 <- fit_shift(D$y, D$phy,ignore.cov=F)   # omit shifts arg for constant rate
fit1 <- fit_shift(D$y, D$phy, shifts, ignore.cov=F)
fit2 <- fit_censored(D$y, D$phy, shifts, ignore.cov=F)

#-------------------------------# 
D <- dataset_elongation(fulltree)
shifts <- find_clades(D$phy, fulltree, c("Serpentes", "Macrostomata"))
bt <- sort(ape::branching.times(D$phy), decreasing=T)
dd <- getDescendantMatrix(D$phy)

bt[as.character(shifts)]

isShiftConfigValid(dmat, bt, shifts)






