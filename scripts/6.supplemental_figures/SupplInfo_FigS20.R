setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

library(ape)
library(cutphylo)
library(scales)

alldat <- read.csv('./data/alldat.csv')

tr <- read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
tr <- ladderize(tr, right = TRUE)

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
tr <- rotate(tr, nn)

tr <- read.tree(text = write.tree(tr)) # otherwise, node rotation doesn't change structure


snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']
snakes <- intersect(snakes, tr$tip.label)
lizards <- intersect(lizards, tr$tip.label)


major_clades = function(phy)
{
    # phy should be the fulltree
    D = alldat
    cl = grep("clade_", colnames(D), value=TRUE)
    cl = cl[-match("clade_tropicalDipsadines", cl)]
    nn = setNames(integer(length(cl)), 
        sapply(strsplit(cl, "clade_", fixed=TRUE), "[[", 2))

    for (i in 1:length(cl))
    {
        tips = intersect(phy$tip.label, D$treename[ D[[ cl[i] ]] == 1L ])
		
		if (length(tips) > 0) {
	        nn[i] = ape::getMRCA(phy, tips)
    	} else {
    		nn[i] <- NA
    	}    
    }
	
	nn <- nn[!is.na(nn)]
    nn
}


spanning_pair_label2 = function(mrca, phy, major_clades)
{
    p = ape::nodepath(phy, mrca, ape::Ntip(phy)+1L)
    cl = intersect(p, major_clades)[1]
    cl = names(major_clades)[match(cl, major_clades)]
    if (!mrca %in% major_clades) {
    	cl <- paste0('within ', cl)
    }
    cl <- paste0(cl, ' (', round(ape::branching.times(phy)[as.character(mrca)], 1), ' mya)')
	
	cl
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

iguanidCol <- '#1b9e77'

# ###########################
# All diet categories

pdf('~/Downloads/figS20.pdf', width = 9, height = 9)
# par(mfrow = c(2,2))

layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE), widths = c(0.55, 0.45))

dat <- read.csv('./data/3.trait_preparation_diversification/dietModeling/diet-proportions-phylo.csv', row.names = 1)
dat <- data.matrix(dat)
dat <- sweep(log(dat), 1, rowMeans(log(dat))) # clr


common <- intersect(rownames(dat), tr$tip.label)
dat <- dat[common,]
dat <- as.matrix(dat)
tree2 <- keep.tip(tr, common)

# partition sum of squares around the arithmetic mean

K <- 25

zz = cut(tree2, dat, K)

n <- 20

s1 <- summary(zz, phy = tree2)
shifts <- s1[, 'partition']
R2 <- s1[, 'R2']

mc = major_clades(tree2)
l <- sapply(shifts[2:(n+1)], function(x) spanning_pair_label2(x, tree2, mc))
    
# do nodes belong to lizards or snakes?
snakeMRCA <- getMRCA(tree2, intersect(tree2$tip.label, snakes))
allSnakeNodes <- c(major_clades(tree2)['Serpentes'], phangorn::Descendants(tree2, major_clades(tree2)['Serpentes'], type = 'all'))
allIguanidNodes <- c(major_clades(tree2)['Iguania'], phangorn::Descendants(tree2, major_clades(tree2)['Iguania'], type = 'all'))
shiftCols <- rep('slateblue1', length(shifts))
shiftCols[shifts %in% allSnakeNodes] <- 'coral2'
shiftCols[shifts %in% allIguanidNodes] <- iguanidCol
table(shiftCols)


par(mar = c(4,4,3,1))

plot(x = 1:(n+1), y = R2[1:(n+1)], type = 'b', pch = rep('', (n+1)), las = 1, bty = 'l', col = 'black', ylim = c(0,1), xpd = NA, xlab = 'partitions', ylab = 'variance explained')

text(1:(n+1), R2[1:(n+1)], c('', as.character(1:n)), col = shiftCols[1:(n+1)], font = 2)

# if (R2[2] < 0.6) {
	# legend("topleft", legend=l, bty="n")
# } else {
	# legend("bottomright", legend=l, bty="n")
# }

# title('all squamates', adj = 0)


edgeCols <- rep('slateblue1', Nedge(tree2))
edgeCols[tree2$edge[,2] %in% phangorn::Descendants(tree2, major_clades(tree2)['Serpentes'], type = 'all')] <- 'coral2'
edgeCols[tree2$edge[,2] %in% phangorn::Descendants(tree2, major_clades(tree2)['Iguania'], type = 'all')] <- iguanidCol
table(edgeCols)

par(mar = c(0,0,0,0))

plot(tree2, edge.color = edgeCols, show.tip.label = FALSE, edge.width = 0.75)
nodelabels(node = shifts[2:(n+1)], pch = 21, bg = 'white', cex = 1.9)

text(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[shifts[2:(n+1)]],
	get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[shifts[2:(n+1)]],
	labels = 1:n, cex = 0.6, font = 2)

legend(x = 10, y = Ntip(tree2)/2, legend = c('lizards', 'iguanian lizards', 'snakes'), pt.bg = c('slateblue1', iguanidCol, 'coral2'), pch = 21, bty = 'n')



## Lizards only

dat <- read.csv('./data/3.trait_preparation_diversification/dietModeling/diet-proportions-phylo.csv', row.names = 1)

dat <- data.matrix(dat)
dat <- sweep(log(dat), 1, rowMeans(log(dat))) # clr

common <- Reduce(intersect, list(rownames(dat), tr$tip.label, lizards))
dat <- dat[common,]
dat <- as.matrix(dat)
tree3 <- keep.tip(tr, common)

K <- 10

zz = cut(tree3, dat, K)

s1 <- summary(zz, phy = tree3)
shifts <- s1[, 'partition']
# shifts <- shifts[-1]
R2 <- s1[, 'R2']

n <- 5


mc = major_clades(tree3)
l <- sapply(shifts[2:(n+1)], function(x) spanning_pair_label2(x, tree3, mc))

# convert partition nodes to nodes in tree2 (all diet data)
shifts2 <- rep(NA, length(shifts))
for (i in 1:length(shifts)) {
	if (!is.na(shifts[i])) {
		shifts2[i] <- getMRCA(tree2, geiger::tips(tree3, shifts[i]))
	}
}
    
# do nodes belong to lizards or snakes?
allIguanidNodes <- c(major_clades(tree2)['Iguania'], phangorn::Descendants(tree2, major_clades(tree2)['Iguania'], type = 'all'))
shiftCols <- rep('slateblue1', length(shifts))
shiftCols[shifts2 %in% allSnakeNodes] <- 'coral2'
shiftCols[shifts2 %in% allIguanidNodes] <- iguanidCol
table(shiftCols)



par(mar = c(4,4,3,1))

plot(x = 1:(n+1), y = R2[1:(n+1)], type = 'b', pch = rep('', (n+1)), las = 1, bty = 'l', col = 'black', ylim = c(0,1), xpd = NA, xlab = 'partitions', ylab = 'variance explained')

text(1:(n+1), R2[1:(n+1)], c('', as.character(1:n)), col = shiftCols, font = 2)



# if (R2[2] < 0.6) {
	# legend("topleft", legend=l, bty="n")
# } else {
	# legend("bottomright", legend=l, bty="n")
# }

# title('lizards only', adj=0)



edgeCols <- rep('slateblue1', Nedge(tree2))
edgeCols[tree2$edge[,2] %in% phangorn::Descendants(tree2, major_clades(tree2)['Serpentes'], type = 'all')] <- 'gray90'
edgeCols[tree2$edge[,2] %in% phangorn::Descendants(tree2, major_clades(tree2)['Iguania'], type = 'all')] <- iguanidCol
table(edgeCols)

par(mar = c(0,0,0,0))

plot(tree2, edge.color = edgeCols, show.tip.label = FALSE, edge.width = 0.75)
nodelabels(node = shifts2[2:(n+1)], pch = 21, bg = 'white', cex = 1.9)

text(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[shifts2[2:(n+1)]],
	get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[shifts2[2:(n+1)]],
	labels = 1:n, cex = 0.6, font = 2)

legend(x = 10, y = Ntip(tree2)/2, legend = c('lizards', 'iguanian lizards'), pt.bg = c('slateblue1', iguanidCol), pch = 21, bty = 'n')


text(grconvertX(0.03, from = 'ndc'), grconvertY(0.97, from = 'ndc'), 'A.', font = 2, cex = 1.3, xpd = NA)
text(grconvertX(0.53, from = 'ndc'), grconvertY(0.97, from = 'ndc'), 'B.', font = 2, cex = 1.3, xpd = NA)
text(grconvertX(0.03, from = 'ndc'), grconvertY(0.465, from = 'ndc'), 'C.', font = 2, cex = 1.3, xpd = NA)
text(grconvertX(0.53, from = 'ndc'), grconvertY(0.465, from = 'ndc'), 'D.', font = 2, cex = 1.3, xpd = NA)

text(grconvertX(0.15, from = 'ndc'), grconvertY(0.97, from = 'ndc'), 'all squamates', font = 2, cex = 1.3, xpd = NA)
text(grconvertX(0.15, from = 'ndc'), grconvertY(0.465, from = 'ndc'), 'lizards only', font = 2, cex = 1.3, xpd = NA)


dev.off()


