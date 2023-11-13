library(Rphylopars)
library(ape)


rescaleTreeFromASR <- function(phy, dat) {
	
	# use Rphylopars for fast ancestral character reconstruction under BM
	if (anyNA(dat)) {
		stop("contains NA's.")
	}
	
	if (inherits(dat, c('numeric', 'integer'))) {
		# anc <- phytools::fastAnc(phy, dat)
		# edge_start <- anc[as.character(phy$edge[,1])]
		# edge_end <- anc[as.character(phy$edge[,2])]
		# new_bl1 <- abs(edge_end - edge_start)
		
		df <- cbind.data.frame(species = phy$tip.label, trait = dat[phy$tip.label])
		phylo_correlated <- TRUE # for some reason, this needs to be set for the univariate to work
	} else {
		df <- cbind.data.frame(species = phy$tip.label, dat[phy$tip.label,])
		phylo_correlated <- FALSE
		
	}
	
	anc <- phylopars(trait_data = df, tree = phy, model = 'BM', phylo_correlated = phylo_correlated)
	head(anc$anc_recon, 2)
	
	# replace tip names with tip node numbers
	## Not currently needed
	if (!identical(phy$tip.label, rownames(anc$anc_recon)[1:Ntip(phy)])) {
		stop()
		match(rownames(anc$anc_recon)[1:Ntip(phy)], phy$tip.label)
	}
	rownames(anc$anc_recon)[1:Ntip(phy)] <- as.character(1:Ntip(phy))
	
	new_bl <- numeric(nrow(phy$edge))
	for (i in 1:nrow(phy$edge)) {
		new_bl[i] <- as.numeric(abs(dist(anc$anc_recon[as.character(phy$edge[i,]),])))
	}
	
	new_bl[is.na(new_bl)] <- 0
	tree_rescaled <- phy
	tree_rescaled$edge.length <- as.numeric(new_bl)
	
	# ## version with mvMORPH::mvgls	
	# fit <- mvMORPH::mvgls(Y ~ 1, data = list(Y = as.matrix(dat[phy$tip.label,])), tree = phy, model = 'BM')
	# anc <- mvMORPH::ancestral(fit)
	
	# # add tip states
	# tipstates <- as.matrix(dat[phy$tip.label,])
	# tipstates <- tipstates[match(phy$tip.label, rownames(tipstates)),]
	# rownames(tipstates) <- paste0('node_', 1:nrow(tipstates))
	# anc <- rbind(tipstates, anc)
	
	# new_bl <- numeric(nrow(phy$edge))
	# for (i in 1:nrow(phy$edge)) {
		# new_bl[i] <- as.numeric(abs(dist(anc[paste0('node_', phy$edge[i,]),])))
	# }
	
	# calculate total root-to-tip distance
	spEdges <- ape::nodepath(tree_rescaled)
	names(spEdges) <- tree_rescaled$tip.label
	spEdges <- lapply(spEdges, function(x) setdiff(match(x, tree_rescaled$edge[,2]), NA))
	rootTipDists <- sapply(spEdges, function(x) sum(tree_rescaled$edge.length[x]))

	
	return(list(tree_rescaled, rootTipDists))
	
}


colorBranches <- function(tree, nodes, nodeCols = 'coral2', defaultCol = gray(0.65)) {

	edge.color <- rep(defaultCol, nrow(tree$edge))
	for (i in 1:length(nodes)) {
		edge.color[which(tree$edge[,1] %in% c(nodes[i], phangorn::Descendants(tree, nodes[i], type = 'all')))] <- nodeCols[i]
	}	

	return(edge.color)
}




setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

alldat <- read.csv('./data/alldat.csv')

tree <- read.tree('./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre')
tree <- read.tree(text = write.tree(ladderize(tree)))

snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']


## DIET
# ancDistDiet

x = read.csv('./data/3.trait_preparation_diversification/dietModeling/diet-proportions-phylo.csv', row.names=1)
x = data.matrix(x)

commonsp = intersect(tree$tip.label, rownames(x))

subtree = keep.tip(tree, commonsp)
x <- x[subtree$tip.label, ]
Y = sweep(log(x), 1, rowMeans(log(x)))

dietASR <- rescaleTreeFromASR(subtree, Y)

# ----------------------
## ELONGATION
elong <- setNames(alldat$elongationIndex, alldat$treename)
elong <- elong[!is.na(elong)]

commonsp <- intersect(tree$tip.label, names(elong))

subtree <- keep.tip(tree, commonsp)
elong <- elong[subtree$tip.label]

elongASR <- rescaleTreeFromASR(subtree, elong)

# ----------------------
## SKULL
x = read.csv("./data/3.trait_preparation_diversification/2D_Adult_Skull_Coord_treenames.csv")
rownames(x) <- x$treename
goodcol = grep('Proc', colnames(x), value = TRUE, ignore.case = TRUE)

commonsp = intersect(tree$tip.label, rownames(x))

subtree = keep.tip(tree, commonsp)
x = x[commonsp, goodcol]

skullASR <- rescaleTreeFromASR(subtree, x)


# ----------------------
## VERTEBRAL
vert <- setNames(alldat$numberPresacralVert, alldat$treename)
vert <- vert[!is.na(vert)]
vert <- log(vert)

commonsp <- intersect(tree$tip.label, names(vert))

subtree <- keep.tip(tree, commonsp)
vert <- vert[subtree$tip.label]

vertASR <- rescaleTreeFromASR(subtree, vert)






pdf('~/Downloads/figS7.pdf', width = 8, height = 8)
par(mfrow = c(2,2))

edgeWidth <- 0.2
letterHeight <- 1.8

snakeNode <- getMRCA(dietASR[[1]], intersect(dietASR[[1]]$tip.label, snakes))
branchCols <- colorBranches(dietASR[[1]], nodes = snakeNode, nodeCols = 'coral2', defaultCol = 'slateblue1')
plot.phylo(dietASR[[1]], show.tip.label = FALSE, type = 'fan', open.angle = 10, edge.width = edgeWidth, edge.color = branchCols, lend = 1, no.margin = TRUE)
mtext('multivariate diet', side = 3, line = -1.5, font = 2)
mtext('A', adj = 0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)


legend('topright', legend = c('lizards', 'snakes'), fill = c('slateblue1', 'coral2'), border = NA, bty = 'n', cex = 1.3)

snakeNode <- getMRCA(elongASR[[1]], intersect(elongASR[[1]]$tip.label, snakes))
branchCols <- colorBranches(elongASR[[1]], nodes = snakeNode, nodeCols = 'coral2', defaultCol = 'slateblue1')
plot.phylo(elongASR[[1]], show.tip.label = FALSE, type = 'fan', open.angle = 10, edge.width = edgeWidth, edge.color = branchCols, lend = 1, no.margin = TRUE)
mtext('elongation', side = 3, line = -1.5, font = 2)
mtext('B', adj = 0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

snakeNode <- getMRCA(skullASR[[1]], intersect(skullASR[[1]]$tip.label, snakes))
branchCols <- colorBranches(skullASR[[1]], nodes = snakeNode, nodeCols = 'coral2', defaultCol = 'slateblue1')
plot.phylo(skullASR[[1]], show.tip.label = FALSE, type = 'fan', open.angle = 10, edge.width = edgeWidth, edge.color = branchCols, lend = 1, no.margin = TRUE)
mtext('multivariate skull shape', side = 3, line = -1.5, font = 2)
mtext('C', adj = 0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

snakeNode <- getMRCA(vertASR[[1]], intersect(vertASR[[1]]$tip.label, snakes))
branchCols <- colorBranches(vertASR[[1]], nodes = snakeNode, nodeCols = 'coral2', defaultCol = 'slateblue1')
plot.phylo(vertASR[[1]], show.tip.label = FALSE, type = 'fan', open.angle = 10, edge.width = edgeWidth, edge.color = branchCols, lend = 1, no.margin = TRUE)
mtext('presacral vertebral count', side = 3, line = -1.5, font = 2)
mtext('D', adj = 0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

dev.off()


