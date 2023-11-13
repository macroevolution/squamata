basedir <- '~/Dropbox/Oz_Crown_Ages/dataArchive'


library(cutphylo) # [remotes|devtools]::install_github("blueraleigh/cutphylo")


ALLDATA = function()
{
	
    BASEDIR = basedir
    TRAITDIR = "scripts/data"
    TRAITFILE = "alldat.csv"
    PATH = file.path(BASEDIR,TRAITDIR,TRAITFILE)
	D <- read.csv(PATH)
	D
	
}

dataset_skull = function(fulltree)
{
    BASEDIR = basedir
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
    return (list(y=y, phy=phy))
}

dataset_diet = function(fulltree)
{
    BASEDIR = basedir
    TRAITDIR = "./data/3.trait_preparation_diversification/dietModeling"
    TRAITFILE = "diet-proportions-phylo.csv"
    PATH = file.path(BASEDIR,TRAITDIR,TRAITFILE)
    x = data.matrix(read.csv(PATH, row.names=1))
    y = sweep(log(x), 1, rowMeans(log(x)))
    common = intersect(rownames(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = y[phy$tip.label, ]
    return (list(y=y, phy=phy))
}

dataset_vertebral = function(fulltree)
{
    D = ALLDATA()
    y = setNames(log(D[[ "numberPresacralVert" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_mass = function(fulltree)
{
    D = ALLDATA()
    y = setNames(log(D[[ "mass" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_elongation = function(fulltree)
{
    D = ALLDATA()
    y = setNames(log(D[[ "elongationIndex" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_climate = function(fulltree)
{
    D = ALLDATA()
    y = data.matrix(log(D[, c("bio7", "bio12", "npp")]))
    rownames(y) = D[[ "treename" ]]
    y = y[apply(y, 1, Negate(anyNA)), ]
    common = intersect(rownames(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = y[phy$tip.label, ]
    return (list(y=y, phy=phy))
}

dataset_bamm = function(fulltree)
{
    D = ALLDATA()
    y = setNames(log(D[[ "bamm" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_clads = function(fulltree)
{
    D = ALLDATA()
    y = setNames(log(D[[ "meanImputedCLADS" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_dietancdist = function(fulltree)
{
    D = ALLDATA()
    y = setNames(D[[ "ancDistDiet" ]], D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_skullancdist = function(fulltree)
{
    D = ALLDATA()
    y = setNames(D[[ "skullAncDist" ]], D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_dietrate = function(fulltree)
{
    D = ALLDATA()
    y = setNames(log(D[[ "dietRate" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

dataset_skullrate = function(fulltree)
{
    D = ALLDATA()
    y = setNames(log(D[[ "skullRate" ]]), D[[ "treename" ]])
    y = y[!is.na(y)]
    common = intersect(names(y), fulltree$tip.label)
    phy = ape::keep.tip(fulltree, common)
    y = t(t(y[phy$tip.label]))
    return (list(y=y, phy=phy))
}

spanning_pair = function(node, phy)
{
    if (node <= ape::Ntip(phy))
    {
        tip = phy$tip.label[node]
        return (c(tip, tip))
    }
    tips = ape::extract.clade(phy, node)$tip.label
    c(tips[1], tail(tips, 1))
}

summary.cutphylo = function(x, ...)
{
    k = length(x$RSS)
    ans = vector("list", k)
    args = list(...)
    phy = args$phy
    stopifnot(!is.null(phy))
    lsp = character(k)
    rsp = character(k)
    for (i in 1:k)
    {
        tmp = spanning_pair(x$cuts[i], phy)
        lsp[i] = tmp[1]
        rsp[i] = tmp[2]
    }
    SST = x$RSS[1]
    RSS = x$RSS
    ESS = SST - RSS
    R2 = signif(ESS / SST, 3)
    P = x$cuts
    P[1] = NA_integer_
    data.frame(
        partition=P
        , left=lsp
        , right=rsp
        , RSS=RSS
        , ESS=ESS
        , R2=R2
    )
}

setwd(basedir)

fulltree = ape::read.tree("./dataArchive/data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")

DATASET = c(
    "clads"
    , "bamm"
    , "diet"
    , "skull"
    , "mass"
    , "vertebral"
    , "elongation"
    , "climate"
    , "dietrate"
    , "skullrate"
    , "dietancdist"
    , "skullancdist"
)

K = 15

outdir <- paste0(basedir, '/dataArchive/data/4.comparative_analyses')

for (d in DATASET)
{
    D = get(sprintf("dataset_%s", d))(fulltree)
    z1 = cut(D$phy, D$y, K, option="arithmetic")
    # z2 = cut(D$phy, D$y, K, option="phylogenetic")
    s1 = summary(z1, phy=D$phy)
    # s2 = summary(z2, phy=D$phy)
    write.csv(s1, file=sprintf(paste0(outdir, '/a-partition/a-partition-R2-%s.csv'), d))
}
