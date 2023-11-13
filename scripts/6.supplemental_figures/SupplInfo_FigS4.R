
## Plot tree with fossil-calibrated nodes indicated

require(MCMCtreeR)

require(phyloch) # stratigraphic legend
data(strat2012)

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

alldat <- read.csv('./data/alldat.csv')

fossilfile <- './data/2.time_calibration_imputation/fossildat.csv'
mcmcPostFile <- './data/2.time_calibration_imputation/MCMCTreeAnalyses/mcmctree_useData2/iterA/mcmc.txt'
mcmctreefile <- './data/2.time_calibration_imputation/MCMCTreeAnalyses/mcmctree_useData2/iterA/FigTree.tre'

mcmcPost <- data.table::fread(mcmcPostFile, data.table = FALSE)

tr1 <- readMCMCtree(mcmctreefile)


timeTree <- tr1[[1]]
timeTree <- ladderize(timeTree)
timeTree <- read.tree(text = write.tree(timeTree, file = ''))

nodeTimes <- cbind.data.frame(oldNodes = rownames(tr1[[2]]), newNodes = NA, tr1[[2]])
nodeTimes$newNodes <- sapply(rownames(tr1[[2]]), function(x) getMRCA(timeTree, geiger::tips(tr1[[1]], x)))

# for each node in tree, calculate the median divergence time from the posterior, and the 95% HPD
for (i in 1:nrow(nodeTimes)) {
	xx <- mcmcPost[, paste0('t_n', nodeTimes[i, 'oldNodes'])]
	nodeTimes[i, 'median'] <- median(xx)
	nodeTimes[i, 'hpd5'] <- HDInterval::hdi(xx, credMass = 0.95)['lower']
	nodeTimes[i, 'hpd95'] <- HDInterval::hdi(xx, credMass = 0.95)['upper']
}

# rescale timeTree with median divergence times
newBL <- numeric(length(timeTree$edge.length))
origBL <- numeric(length(timeTree$edge.length))
for (i in 1:nrow(timeTree$edge)) {
	xx <- nodeTimes[nodeTimes$newNodes %in% timeTree$edge[i, ],]
	if (nrow(xx) == 2) {
		origBL[i] <- abs(diff(xx[, 'mean']))
		newBL[i] <- abs(diff(xx[, 'median']))
	} else if (nrow(xx) == 1) {
		origBL[i] <- xx[, 'mean']
		newBL[i] <- xx[, 'median']	
	} else stop()
}

# do we recover the original branch lengths?
table(abs(origBL - timeTree$edge.length) < 0.0000001)
range(origBL - timeTree$edge.length)

timeTree$edge.length <- newBL

timeTree <- geiger:::rescale.phylo(timeTree, model = 'depth', depth = max(branching.times(timeTree)) * 100)


nodeTimes[, 3:8] <- nodeTimes[, 3:8] * 100


# fossil calibration data
fossildat <- read.csv(fossilfile)


# label all families
uniqueFams <- unique(alldat$family)
famNodes <- rep(NA, length(uniqueFams))
names(famNodes) <- uniqueFams
for (i in 1:length(uniqueFams)) {
	
	commonsp <- intersect(alldat[alldat$family == uniqueFams[i], 'treename'], timeTree$tip.label)
	if (length(commonsp) > 1) {
		famNodes[i] <- getMRCA(timeTree, commonsp)
	} else if (length(commonsp) == 1) {
		famNodes[i] <- which(timeTree$tip.label == commonsp)
	}
}
famNodes <- famNodes[!is.na(famNodes)]

maxAge <- max(branching.times(timeTree))


fossildat <- fossildat[which(fossildat$use == 'yes'),]
fossildat <- fossildat[which(fossildat$tree == 'genomicBackbone'),]
# which(sapply(fossildat$spanning.taxa, function(x) all(strsplit(x, split = '\\|')[[1]] %in% timeTree$tip.label)) == TRUE)

# redefine node labels
## and save for use throughout
labelAges <- c()
for (i in 1:nrow(fossildat)) {
	labelAges[i] <- branching.times(timeTree)[as.character(getMRCA(timeTree, strsplit(fossildat[i, 'spanning.taxa'], split = '\\|')[[1]]))]
}
fossildat[order(labelAges, decreasing = TRUE), 'newNodeLabel'] <- paste0('C', 1:nrow(fossildat))



pdf('~/Downloads/figS4.pdf', width = 8, height = 10)

marvec <- c(3,1,1,5)

par(mar = marvec)

# plot.phylo(timeTree, no.margin = F, cex = 0.4, label.offset = 1, edge.color = 'white', tip.color = 'white')
phytools::plotTree(timeTree, fsize = 0.4, mar = marvec, lwd = 0.8, plot = F)

start <- maxAge
end <- maxAge - 10
while (end > 0) {
	rect(xright = start, xleft = end, ytop = Ntip(timeTree) + 1, ybottom = 0, border = NA, col = gray(0.97))
	start <- end - 10
	end <- start - 10
}

phytools::plotTree(timeTree, fsize = 0.4, ftype = 'i', mar = marvec, lwd = 0.9, add = TRUE)
phyloch::axisGeo(GTS = strat2012, unit = c("period"), cex = 0.8)


pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

pp$xx[nodeTimes$newNodes]
pp$yy[nodeTimes$newNodes]

segments(x0 = maxAge - c(50, 100, 150, 200, 250, 300), x1 = maxAge - c(50, 100, 150, 200, 250, 300), y0 = 1, y1 = Ntip(timeTree) + 3, lty = 3, xpd = NA)
text(maxAge - c(50, 100, 150, 200, 250, 300), y = Ntip(timeTree) + 2, labels = c(50, 100, 150, 200, 250, 300), pos = 3, cex = 0.5, xpd = NA)

nn <- sapply(fossildat$spanning.taxa, function(x) getMRCA(timeTree, strsplit(x, '\\|')[[1]]))

nodelabels(node = nn, pch = 21, frame = 'circle', bg = adjustcolor('lightblue', alpha.f = 0.5), cex = 2)

segments(x0 = maxAge - nodeTimes[, 'hpd5'], x1 = maxAge - nodeTimes[, 'hpd95'], y0 = pp$yy[nodeTimes$newNodes], y1 = pp$yy[nodeTimes$newNodes], col = 'dodgerblue2', lwd = 3, lend = 1, xpd = NA)


# nodelabels(text = fossildat$nodeLabel, node = nn, frame = 'circle', bg = adjustcolor('lightblue', alpha.f = 0.9), cex = 0.5)


nodelabels(text = fossildat$newNodeLabel, node = nn, frame = 'none', adj = c(1.5,1.5), cex = 0.5, col = 'darkblue', font = 2)

famNodes <- famNodes[order(pp$yy[famNodes])]

doubleLabels <- intersect(famNodes, nn)
famInd1 <- which(!famNodes %in% doubleLabels)
famInd2 <- which(famNodes %in% doubleLabels)

xbars <- grconvertX(0.87, from = 'ndc')

for (i in 1:length(famNodes)) {
	
	yrange <- range(which(timeTree$tip.label %in% geiger::tips(timeTree, famNodes[i])))
	
	segments(x0 = xbars, x1 = xbars, y0 = yrange[1] - 0.4, y1 = yrange[2] + 0.4, lend = 1, xpd = NA)
	
	text(x = xbars, y = mean(yrange), names(famNodes)[i], pos = 4, cex = 0.45, xpd = NA)
	
}

dev.off()


