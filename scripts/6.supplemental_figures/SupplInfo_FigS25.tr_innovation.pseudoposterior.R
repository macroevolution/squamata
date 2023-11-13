setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

library(phylo)
library(bm)

source('./scripts/3.trait_preparation_diversification/calcNetInnovation.R')

alldat <- read.csv('./data/alldat.csv')

# location of pseudoposterior of trees
trees = ape::read.tree("./data/1.tree_inference/pseudoposterior/pseudoposterior.100.trees")

# Function to calculate TR
calcTR <- function(dat, phy, log = TRUE) {
	
	if (!is.matrix(dat)) {
		dat <- data.matrix(dat)	
	}
	
	dat <- dat[complete.cases(dat),, drop = FALSE]
	
	commonTaxa <- intersect(rownames(dat), phylo::tiplabels(phy))
	traitTree <- keep.tip(phy, commonTaxa)
	dat <- dat[tiplabels(traitTree),, drop = FALSE]
	
	if (log) {
		dat <- log(dat[,, drop = FALSE])
	}
		

	# univariate
	if (ncol(dat) == 1) {
		
		datRate <- bm.tiprate(setNames(dat[,1], rownames(dat)), traitTree)
		names(datRate) <- rownames(dat)
		
	} else {
	
		# multivariate
		rates <- bm.mvtiprate(dat, traitTree)

		datRate <- numeric(nrow(dat))
		names(datRate) <- rownames(dat)
		for (ii in 1:length(rates)) {
			datRate[ii] <- sum(diag(rates[[ii]]))
			# eig <- eigen(rates[[ii]])
			# climateNicheRate[ii] <- sqrt(max(eig$values))
		}
	}
	
	return(datRate)
	
}

bdir = '~/Dropbox/Oz_Crown_Ages/dataArchive'
source("./scripts/4.comparative_analyses/model-fits/fxns-datasets.R")

res = vector("list", length(trees) * 4)
for (j in 1:length(trees)) {
  
    treenum = j
		
		phy = trees[[j]]
		# mike's script assumes a unique format
		phy_mg = phylo::read.newick(text = ape::write.tree(phy))
		DATASET = list(
		  dataset_skull(phy, BASEDIR = bdir),
		  dataset_diet_noTree(phy, BASEDIR = bdir),
		  dataset_mass(phy, BASEDIR = bdir),
		  dataset_elongation(phy, BASEDIR = bdir)
		)
		names(DATASET) = c("skull", "diet_noTree",
		                   "mass", "elongation")
		# LOG FALSE, LOG TRUE
		# LOG TRUE, LOG TRUE
		
		mass = as.vector(DATASET[["mass"]]$y)
		names(mass) = rownames(DATASET[["mass"]]$y)
		DATASET[["mass"]]$y = mass
		
		elong = as.vector(DATASET[["elongation"]]$y)
		names(elong) = rownames(DATASET[["elongation"]]$y)
		DATASET[["elongation"]]$y = elong
	
		for (i in 1:length(DATASET)) {		
		  # keep mahalanobis distances
		  inno = calcNetInnovation(DATASET[[i]]$y, phy_mg)[ ,2]
		  # don't log anything bc it has already been
		  # logged in fxns-datasets.R
		  tr = calcTR(DATASET[[i]]$y, phy_mg, log = FALSE)
		  dd = data.frame(treename = names(inno),
		                  metric = names(DATASET)[[i]],
		                  psi = inno,
		                  tr = tr[names(inno)],
		                  treenum = treenum)
		  rownames(dd) = NULL
		  res[[4 * (j - 1) + i]] = dd
		  cat(i, " * of * ", j, "\n", sep = "")
		}
}
res2 = do.call("rbind", res)
write.csv(res2, "./manuscript/figure_data/tr_innovation.pseudoposterior.csv")
