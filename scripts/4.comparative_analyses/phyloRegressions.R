library(phylo) # remotes::install_github('blueraleigh/phylo')
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


runmodel = function(model, alldat, fulltree, use.brlen=FALSE) 
{    
    if (!use.brlen)
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


BASE_DIR = "~/Dropbox/Oz_Crown_Ages"
setwd(BASE_DIR)

alldat <- read.csv('./dataArchive/data/alldat.csv')


fulltree = ape::read.tree("./dataArchive/data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")

# set certain variables to factor
alldat$prehensionMechanism <- factor(alldat$prehensionMechanism, levels = c('jaw', 'both', 'lingual'))
alldat$parity <- factor(alldat$parity, levels = c('Oviparous', 'Mixed', 'Viviparous'))
alldat$foragingMode <- factor(alldat$foragingMode)
alldat$combinedKinesis <- factor(alldat$combinedKinesis, levels = c('akinetic', 'highKinesis', 'hyperkinesis', 'lowKinesis', 'midKinesis'))
alldat$geogRadiation <- factor(alldat$geogRadiation)

for (i in grep("^geog_", colnames(alldat))) {
	alldat[, i] = factor(alldat[, i])
}

for (i in grep("^clade_", colnames(alldat))) {
	alldat[, i] = factor(alldat[, i])
}



# Model definitions with ClaDS rates
M = list(
  
  M_TROPHIC1 = log(meanImputedCLADS) ~ dietPC1 + dietPC2 + log(dietBreadth)
  
  , M_TROPHIC2 = log(meanImputedCLADS) ~ foragingMode
  
  , M_TROPHIC3 = log(meanImputedCLADS) ~ chemosensory_index 
  
  , M_PARITY = log(meanImputedCLADS) ~ parity
  
  , M_KINESIS = log(meanImputedCLADS) ~ combinedKinesis + prehensionMechanism
  
  , M_BIOCLIM = log(meanImputedCLADS) ~ climPC1 + climPC2
  
  , M_MORPH = log(meanImputedCLADS) ~ log(mass) + log(completeSVL) + log(elongationIndex)
  
  , M_MERISTIC1 = log(meanImputedCLADS) ~ numberDigits + numberLimbs
  
  , M_MERISTIC2 = log(meanImputedCLADS) ~ numberPresacralVert 
  
  , M_SKULL = log(meanImputedCLADS) ~ skullPC1 + skullPC2
  
  , M_GEO = log(meanImputedCLADS) ~ centroidLat + elev
  
)


# Model definitions with BAMM rates
M_bamm = list(
  
  M_TROPHIC1 = log(bamm) ~ dietPC1 + dietPC2 + log(dietBreadth) 
  
  , M_TROPHIC2 = log(bamm) ~ foragingMode
  
  , M_TROPHIC3 = log(bamm) ~ chemosensory_index 
  
  , M_PARITY = log(bamm) ~ parity
  
  , M_KINESIS = log(bamm) ~ combinedKinesis + prehensionMechanism
  
  , M_BIOCLIM = log(bamm) ~ climPC1 + climPC2
  
  , M_MORPH = log(bamm) ~ log(mass) + log(completeSVL) + log(elongationIndex)
  
  , M_MERISTIC1 = log(bamm) ~ numberDigits + numberLimbs
  
  , M_MERISTIC2 = log(bamm) ~ numberPresacralVert 
  
  , M_SKULL = log(bamm) ~ skullPC1 + skullPC2
  
  , M_GEO = log(bamm) ~ centroidLat + elev
  
)

M_FIT = lapply(M, runmodel, alldat, fulltree)

M_FIT_bamm = lapply(M_bamm, runmodel, alldat, fulltree)

lapply(M_FIT, summary)
lapply(M_FIT, function(x) summary(x)$adj.r.squared)
lapply(M_FIT_bamm, summary)
lapply(M_FIT_bamm, function(x) summary(x)$adj.r.squared)

round(cbind(clads = sapply(M_FIT, function(x) summary(x)$adj.r.squared), bamm = sapply(M_FIT_bamm, function(x) summary(x)$adj.r.squared)), 3)

# with branch lengths
M_FIT = lapply(M, runmodel, alldat, fulltree, use.brlen=TRUE)

M_FIT_bamm = lapply(M_bamm, runmodel, alldat, fulltree, use.brlen=TRUE)

