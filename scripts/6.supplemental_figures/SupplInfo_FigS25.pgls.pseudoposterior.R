library(phylo)
#
#
# Defines a bunch of linear models of the form
#
#  Y = XB + e      (1)
#
# Y is a n-by-1 column vector of log-transformed ClaDS tiprates
# X is a n-by-p matrix of predictors
# B is a p-by-1 column vector of predictor coefficients
# e is a n-by-1 column vector of (correlated) residuals
#
# Assumes that the correlation structure of e is
# roughly Brownian. In particular, that the 
# covariance between e[i] and e[j] is proportional
# to the number of shared branches between their
# MRCA and the root of the phylogeny.
#
# We compute a contrast matrix C (n-by-[n minus 1]), and
# premultiply (1) by its transpose, CT, to arrive at a new
# model
#
#  Y' = X'B + e'   (2)
#
#  Y' = CT * Y    [n minus 1]-by-1 column vector of tiprate contrasts
#  X' = CT * X    [n minus 1]-by-p column vector of predictor contrasts
#  e' = CT * e    [n minus 1]-by-1 column vector of i.i.d. residuals
#
# If X included an intercept, premultiplying by CT sets it to zero.
#
# Equation (2) is the form that we use below to estimate B.
#
#
# Binary matrix (ntip-by-nnode) that records membership in a clade
treematrix = function(phy) 
{
  X = matrix(0L, phylo::Ntip(phy), phylo::Nnode(phy))
  for (node in phylo::descendants(phylo::root(phy), phy)) 
  {
    if (node <= phylo::Ntip(phy))
      X[node, node] = 1L
    else
      X[phylo::tips(node, phy), node] = 1L
  }
  X[, phylo::root(phy)] = 1L
  structure(X, dimnames=list(NULL, 1:phylo::Nnode(phy)))
}

# Matrix (ntip-by-[ntip minus 1]) of independent contrast vectors
icmatrix = function(phy) 
{
  X = treematrix(phy)
  C = X
  postorder = phylo::descendants(
    phylo::root(phy), phy, "INTERNAL_NODES_ONLY", "POSTORDER")
  postorder = c(postorder, phylo::root(phy))
  V = phylo::brlens(phy)
  for (p in postorder)
  {
    k = phylo::children(phy, p)
    lf = k[1]
    rt = k[2]
    V[p] = V[p] + (V[lf] * V[rt]) / (V[lf] + V[rt])
    X[, p] = (V[rt] * X[, lf] + V[lf] * X[, rt]) / (V[lf] + V[rt])
    C[, p] =  (X[, lf] - X[, rt]) / sqrt(V[lf] + V[rt]) 
  }
  C[, -(1:phylo::Ntip(phy))]
}



runmodel = function(model, alldat, fulltree) 
{    
  fulltree$edge.length = rep(1, nrow(fulltree$edge))
  phy = phylo::read.newick(text=ape::write.tree(fulltree))
  
  X = model.matrix(model, data=alldat, na.action=na.omit)
  
  response <- all.vars(model)[1]
  
  Y = structure(
    log(alldat[rownames(X), response]), 
    names=alldat[rownames(X), ]$treename)
  
  spp = intersect(phylo::tiplabels(phy), names(Y))
  
  keep = names(Y) %in% spp 
  
  X = X[keep, , drop=FALSE]
  Y = Y[keep]
  
  tree = phylo::keep.tip(phy, spp)
  ord = match(phylo::tiplabels(tree), names(Y))
  
  stopifnot(!anyNA(ord))
  
  X = X[ord, , drop=FALSE]
  Y = Y[ord]
  
  C = icmatrix(tree)
  CT = t(C)
  
  # Premultiplying by transpose of contrast matrix changes
  # our problem from PGLS to OLS (and drastically lowers 
  # computation times)
  #
  #    Stone, E.A. (2011) Why the phylogenetic regression appears robust 
  #    to tree misspecification. Systematic Biology, 60, 245–260.
  #
  # Specifically, equation (7)
  YY = CT %*% Y
  XX = (CT %*% X)[, -1L, drop=FALSE]
  
  colnames(XX) = gsub("^log\\(", "ln", colnames(XX))
  colnames(XX) = gsub("\\)$", "", colnames(XX))
  
  model = as.formula(
    paste0('ln_', response, ' ~ ', paste0(colnames(XX), collapse=" + ")))
  
  model = update(model, . ~ . - 1)
  
  colnames(YY) = paste0('ln_', response)
  
  dat = as.data.frame(cbind(YY, XX))
  
  fit = lm(model, data=dat)
  
  #fit$orig.data = cbind(
  #    data.frame(lnClaDS=Y), 
  #    as.data.frame(X[, -1L, drop=FALSE]))
  
  fit
}

setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

alldat <- read.csv('./data/alldat.csv')

# set certain variables to factor
alldat$prehensionMechanism <- factor(alldat$prehensionMechanism, levels = c('jaw', 'both', 'lingual'))
alldat$parity <- factor(alldat$parity, levels = c('Oviparous', 'Mixed', 'Viviparous'))
alldat$foragingMode <- factor(alldat$foragingMode)
alldat$combinedKinesis <- factor(alldat$combinedKinesis, levels = c('akinetic', 'highKinesis', 'hyperkinesis', 'lowKinesis', 'midKinesis'))

trees = ape::read.tree("./data/1.tree_inference/pseudoposterior/pseudoposterior.100.trees")

# this is: data/2.time_calibration_imputation/pseudoposterior-imputed-trees.zip
tt = list.files("./taxon-imputation/pseudoposterior-full-trees/")

d = read.table("./data/3.trait_preparation_diversification/CLaDS/pseudoposterior_100_imputations.CLaDS_rates.txt", row.names = 1)
colnames(d) = gsub("X", "", colnames(d))

fits = vector("list", length(trees))

M = list(
  # note this is somewhat weird bc diet PC1 & 2 are defined by diet proportions
  # inferred  under a phylogenetic model based on the primary tree
  # - but for the intents here, this will do
  M_TROPHIC1 = log(clads_tmp) ~ dietPC1 + dietPC2 + log(dietBreadth)
  
  , M_TROPHIC2 = log(clads_tmp) ~ foragingMode
  
  , M_TROPHIC3 = log(clads_tmp) ~ chemosensory_index 
  
  , M_PARITY = log(clads_tmp) ~ parity
  
  , M_KINESIS = log(clads_tmp) ~ combinedKinesis + prehensionMechanism
  
  , M_BIOCLIM = log(clads_tmp) ~ climPC1 + climPC2
  
  , M_MORPH = log(clads_tmp) ~ log(mass) + log(completeSVL) + log(elongationIndex)
  
  , M_MERISTIC1 = log(clads_tmp) ~ numberDigits + numberLimbs
  
  , M_MERISTIC2 = log(clads_tmp) ~ numberPresacralVert + numberCaudalVert 
  
  , M_SKULL = log(clads_tmp) ~ skullPC1 + skullPC2
  
  , M_GEO = log(clads_tmp) ~ centroidLat + elev
)

for (i in 1:length(trees)) {
  alldat$clads_tmp = d[ match(alldat$treename, rownames(d)), as.character(i)]
  M_FIT = lapply(M, runmodel, alldat, trees[[i]])
  fits[[i]] = sapply(M_FIT, function(x) summary(x)$adj.r.squared)
  cat(i, "\n")
}
saveRDS(fits, file="pgls.pseudoposterior.Rds")
