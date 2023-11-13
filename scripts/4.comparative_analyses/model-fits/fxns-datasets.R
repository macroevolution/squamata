ALLDATA = function(BASEDIR)
{
	
    TRAITDIR = "data"
    TRAITFILE = "alldat.csv"
    PATH = file.path(BASEDIR,TRAITDIR,TRAITFILE)
	D <- read.csv(PATH)
	D
	
}

dataset_skull = function(fulltree, BASEDIR = file.path(DROPBOX, "Oz_Crown_Ages"))
{
    TRAITDIR = "./data/3.trait_preparation_diversification"
    TRAITFILE = "2D_Adult_Skull_Coord_treenames.csv"
    PATH = file.path(BASEDIR,TRAITDIR,TRAITFILE)
    y = read.csv(PATH)
    rownames(y) <- y$treename
    y <- y[, grep('ProcCoord', colnames(y), value = TRUE)]
    y <- data.matrix(y)
    common = intersect(rownames(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = y[phy$tip.label, ]

    #r = cor(y)
    #diag(r) = 0
    #which(r > .89, arr.ind=TRUE) # 13 and 15, 24 and 40

    return (list(y=y, phy=phy))
}

dataset_diet = function(fulltree, BASEDIR = file.path(DROPBOX, "Oz_Crown_Ages"))
{
    TRAITDIR = "data/3.trait_preparation_diversification/dietModeling"
    TRAITFILE = "diet-proportions-phylo.csv"
    PATH = file.path(BASEDIR,TRAITDIR,TRAITFILE)
    x = data.matrix(read.csv(PATH, row.names=1))
    #y = sweep(log(x), 1, rowMeans(log(x)))
    d = ncol(x)
    y = log(sweep(x[,-d], 1, x[,d], "/"))
    common = intersect(rownames(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = y[phy$tip.label, ]

    #r = cor(y)
    #diag(r) = 0
    #which(r > .89, arr.ind=TRUE) # 6 and 18, 23 and 29

    return (list(y=y, phy=phy))
}

dataset_diet_noTree = function(fulltree, BASEDIR = file.path(DROPBOX, "Oz_Crown_Ages"))
{
    
    TRAITDIR = "data/3.trait_preparation_diversification/dietModeling"
    TRAITFILE = "diet-proportions-nonphylo.csv"
    PATH = file.path(BASEDIR,TRAITDIR,TRAITFILE)
    x = data.matrix(read.csv(PATH, row.names=1))
    #y = sweep(log(x), 1, rowMeans(log(x)))
    d = ncol(x)
    y = log(sweep(x[,-d], 1, x[,d], "/"))
    common = intersect(rownames(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = y[phy$tip.label, ]

    #r = cor(y)
    #diag(r) = 0
    #which(r > .89, arr.ind=TRUE) # 6 and 18, 23 and 29

    return (list(y=y, phy=phy))
}


dataset_dietPC_noTree = function(fulltree, cutoff=0.95)
{
    BASEDIR = file.path(DROPBOX, "Oz_Crown_Ages")
    TRAITDIR = "data/3.trait_preparation_diversification/dietModeling"
    TRAITFILE = "proportions-nonphylo.csv"
    PATH = file.path(BASEDIR,TRAITDIR,TRAITFILE)
    x = data.matrix(read.csv(PATH, row.names=1))
    d = ncol(x)
    y = sweep(log(x), 1, rowMeans(log(x)))
    common = intersect(rownames(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = y[phy$tip.label, ]
    pca = prcomp(y)
    r2 = cumsum(pca$sdev^2/sum(pca$sdev^2))
    k = which(r2 >= cutoff)[1]
    return (list(y=pca$x[,1:k], phy=phy))
}



dataset_vertebral = function(fulltree, BASEDIR)
{
    D = ALLDATA(BASEDIR = BASEDIR)
    y = setNames(log(D[[ "numberPresacralVert" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_mass = function(fulltree, BASEDIR)
{
    D = ALLDATA(BASEDIR = BASEDIR)
    y = setNames(log(D[[ "mass" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_elongation = function(fulltree, BASEDIR)
{
    D = ALLDATA(BASEDIR = BASEDIR)
    y = setNames(log(D[[ "elongationIndex" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_climate = function(fulltree, BASEDIR)
{
    D = ALLDATA(BASEDIR = BASEDIR)
    y = data.matrix(log(D[, c("bio7", "bio12", "npp")]))
    rownames(y) = D[[ "treename" ]]
    y = y[apply(y, 1, Negate(anyNA)), ]
    common = intersect(rownames(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = y[phy$tip.label, ]
    return (list(y=y, phy=phy))
}
