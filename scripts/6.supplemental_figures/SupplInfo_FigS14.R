basedir <- '~/Dropbox/Oz_Crown_Ages/dataArchive'

filedir <- './data/4.comparative_analyses/a-partition'

ALLDATA = function()
{
	
    BASEDIR = basedir
    TRAITDIR = "./data"
    TRAITFILE = "alldat.csv"
    PATH = file.path(BASEDIR,TRAITDIR,TRAITFILE)
	D <- read.csv(PATH)
	D
	
}

major_clades = function(phy)
{
    # phy should be the fulltree
    D = ALLDATA()
    cl = grep("clade_", colnames(D), value=TRUE)
    cl = cl[-match("clade_tropicalDipsadines", cl)]
    nn = setNames(integer(length(cl)), 
        sapply(strsplit(cl, "clade_", fixed=TRUE), "[[", 2))

    for (i in 1:length(cl))
    {
        tips = intersect(phy$tip.label, D$treename[ D[[ cl[i] ]] == 1L ])

        nn[i] = ape::getMRCA(phy, tips)
        
    }

    nn
}

spanning_pair_label2 = function(pair, phy, major_clades)
{
    p = ape::nodepath(fulltree, ape::getMRCA(phy, pair), ape::Ntip(phy)+1L)
    cl = intersect(p, major_clades)[1]
    cl = names(major_clades)[match(cl, major_clades)]
    if (!ape::getMRCA(phy, pair) %in% major_clades) {
    	cl <- paste0('within ', cl)
    }
    cl <- paste0(cl, ' (', round(ape::branching.times(phy)[as.character(ape::getMRCA(phy, pair))], 1), ' mya)')
	
	cl
}


foo = function(p, n, dataset, plot_args, ...)
{
    shifts = p$partition[-1]
    text = as.character(seq_along(shifts))

    mc = major_clades(fulltree)
    l = unname(apply(p[2:(n+1), ], 1, function(r) {
        spanning_pair_label2(r[2:3], fulltree, mc)
    }))
    
    # do nodes belong to lizards or snakes?
    D = ALLDATA()
    allsnakes <- D[D$clade_Serpentes == 1, 'treename']
	shiftCols <- rep('slateblue1', length(shifts))
	for (i in 1:nrow(p)) {
		if (all(p[i, c('left', 'right')] %in% allsnakes)) {
			shiftCols[i] <- 'coral2'
		}
	}

    for (i in 1:n)
        l[i] = sprintf("%d. %s", i, l[i])
 
    plot_args = c(list(
        1:(n+1)
        , p$R2[1:(n+1)]
        , type='b'
        , pch=rep('', (n+1))
        , las=1
        , bty="l"
        , col = 'black'
        , ylim=c(0,1)), plot_args)
        
    do.call(plot, plot_args)
    
    text(1:(n+1), p$R2[1:(n+1)], c('', as.character(1:n)), col = shiftCols[1:(n+1)], font = 2)
    

    if (p$R2[2] < 0.6)
        legend("topleft", legend=l, bty="n", ...)
    else
        legend("bottomright", legend=l, bty="n", ...)



    title(dataset, adj=0)
}

TREEFILE = file.path(basedir, "./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")

fulltree = ape::read.tree(TREEFILE)

DATASET = c(
    "clads"
    , "bamm"
    , "dietrate"
    , "skullrate"
    , "dietancdist"
    , "skullancdist"
    , "diet"
    , "skull"
    , "mass"
    , "vertebral"
    , "elongation"
    , "climate"
    
)


DATASET_TITLE = setNames(c(
    "A. Speciation rate (CLaDS)"
    , "B. Speciation rate (BAMM)"
    , bquote(bold("C."~"TR"["DIET"]))
    , bquote(bold("D."~TR["SKULL"]))
    , bquote(bold("E."~psi["DIET"]))
    , bquote(bold("F."~psi["SKULL"]))
    , "G. multivariate diet"
    , "H. multivariate skull shape"
    , "I. Mass"
    , "J. Presacral vertebral count"
    , "K. Elongation index"
    , "L. Multivariate climate"
), DATASET)

K = 5

setwd(basedir)

pdf(file="~/Downloads/figS14.pdf", width=8.5, height=11)
par(mfrow=c(4, 3), xpd=NA, mar=c(2,2,2,2), oma=c(3,3,3,3))
for (d in DATASET)
{
    p = read.csv(sprintf(paste0(filedir, '/a-partition-R2-%s.csv'), d), row.names=1)
    foo(p, K, DATASET_TITLE[[d]], plot_args=list(xlab="", ylab=""), cex=0.75)
    i = match(d, DATASET)
    if (i %in% c(1,4,7,10))
        mtext("Variance explained", 2, 2.5, cex=0.8)
    if (i %in% c(10,11,12))
        mtext("n partitions", 1, 2.5, cex=0.8)
}
dev.off()

