# Demonstration of how we calculated the phenotypic tip rate (TR) for univariate and multivariate traits
# Demonstration of how we calculated net innovation

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

library(phylo) # remotes::install_github('blueraleigh/phylo')
library(bm) # remotes::install_github('blueraleigh/bm')

treefile <- './data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre'
traitfile <- './data/alldat.csv'
skullShapeFile <- './data/3.trait_preparation_diversification/2D_Adult_Skull_Coord_treenames.csv'

tree <- read.newick(treefile)
alldat <- read.csv(traitfile)

###################################
# TR for log univariate SVL

svl <- setNames(alldat$completeSVL, alldat$treename)
svl <- svl[!is.na(svl)]

commonTaxa <- intersect(names(svl), tiplabels(tree))
traitTree <- keep.tip(tree, commonTaxa)
svl <- svl[tiplabels(traitTree)]
svlRate <- bm.tiprate(log(svl), traitTree)
names(svlRate) <- names(svl)

head(svlRate)


##################################
# TR for multivariate skull shape

skull <- read.csv(skullShapeFile)
rownames(skull) <- skull$treename
skull <- skull[, grep('ProcCoord', colnames(skull), value = TRUE)]

commonTaxa <- intersect(rownames(skull), tiplabels(tree))
traitTree <- keep.tip(tree, commonTaxa)
skull <- skull[tiplabels(traitTree), ]

skullRate <- bm.mvtiprate(skull, traitTree)
skullRate <- sapply(skullRate, function(x) sum(diag(x)))
names(skullRate) <- rownames(skull)




####################################
# Net Innovation

source('./scripts/3.trait_preparation_diversification/calcNetInnovation.R')

## this function internally subsets to taxa common with the tree
## two distances are returned: euclidian and mahalanobis (we use mahalanobis as it is standardized and therefore comparable across traits)

# net innovation for SVL

svl <- setNames(alldat$completeSVL, alldat$treename)
svl <- svl[!is.na(svl)]

ancdistSVL <- calcNetInnovation(log(svl), tree)

head(ancdistSVL)


# for a multivariate dataset like skull shape

skull <- read.csv(skullShapeFile)
rownames(skull) <- skull$treename
skull <- skull[, grep('ProcCoord', colnames(skull), value = TRUE)]

ancdistSkull <- calcNetInnovation(skull, tree)

head(ancdistSkull)






