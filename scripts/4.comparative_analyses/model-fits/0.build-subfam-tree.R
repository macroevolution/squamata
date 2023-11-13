library(ape)
 
# From base path:

shift_path <- "/Users/drabosky/Dropbox (University of Michigan)/manuscriptwork/SQCL/bigphylo/Oz_Crown_Ages/scripts/fit-shifts"

setwd("/Users/drabosky/Dropbox (University of Michigan)/manuscriptwork/SQCL/bigphylo/Oz_Crown_Ages")
source("scripts/empiricalScripts/6.generateFullDataset.R")
setwd(shift_path)

source("fxns-find-clades.R")
 
# alldat is the master data 
#------------------------------------#

DROPBOX <- "/Users/drabosky/Dropbox (University of Michigan)/manuscriptwork/SQCL/bigphylo"
BASEDIR <- file.path(DROPBOX, "Oz_Crown_Ages")
 
TREEFILE <- file.path(BASEDIR, "./final-trees/best_ultrametric_fulltree_ddBD_revision.tre")
v <- ape::read.tree(TREEFILE)

#------------------------------------# 
  
# Pare full tree down to subfamily-level clades. This is the backbone
#   that defines the set of possible shift nodes

sum(is.na(alldat$subfamily))
xx <- alldat[!is.na(alldat$subfamily), ]

#-------------------------------#
# Any subfamilies missing from tree?

in_both <- intersect(xx$treename, v$tip.label)
v       <- ape::drop.tip(v, setdiff(v$tip.label, in_both))
xx      <- xx[xx$treename %in% in_both, ]
  
 
#-------------#
uf <- unique(xx$subfamily)

fspan <- matrix("", nrow= length(uf), ncol=2)
rownames(fspan) <- uf

tipcounts <- matrix(NA, nrow=length(uf), ncol=2)

for (ii in 1:length(uf)){
	
	iset <- intersect(v$tip.label, xx$treename[xx$subfamily == uf[ii]])
	tipcounts[ii, 1] <- length(iset)
	if (length(iset) >= 2){
		xc <- extract.clade(v, node = getMRCA(v, iset))
		xc <- drop.tip(xc, setdiff(xc$tip.label, iset))
		fspan[ii, ] <- getSpanningTips(xc)
		tipcounts[ii, 2] <- length(xc$tip.label)
	}
	
	if (length(iset) == 1){
		fspan[ii, 1] <- iset
		tipcounts[ii, 2] <- 1
	}
	
}

tmp <- as.vector(fspan)
tmp <- tmp[tmp != ""]
vsf <- drop.tip(v, setdiff(v$tip.label, tmp))

cmat <- buildSpanningSet_DF(vsf)

write.tree(vsf, file = "subfamily-pairs-squams-r1.tre")

#----------------------------------------------#
# subfamily-pairs-squams.tre is the base tree 
 


