# speciation rates figure

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

vertical <- TRUE

require(ape)
require(classInt)
require(phytools)
require(scales)

source('./scripts/6.supplemental_figures/violin.R')

alldat <- read.csv('./data/alldat.csv')

nodeRotationOption <- 1 # good for vertical
# nodeRotationOption <- 2 # good for horizontal

tr <- read.tree('./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre')
tr <- ladderize(tr, right = ifelse(nodeRotationOption == 1, TRUE, FALSE))

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
tr <- rotate(tr, nn)

## This one does NOT need to be rotated if ladderizing with right = FALSE!
# rotate node to avoid high skull innovation chameleons being right next to snakes.
if (nodeRotationOption == 2) {
	tr <- rotate(tr, getMRCA(tr, grep('Uromastyx|Brookesia', tr$tip.label, value = TRUE)))
}

if (nodeRotationOption == 2) {
	tr <- rotate(tr, getMRCA(tr, grep('Saara|Draco', tr$tip.label, value = TRUE)))
}


tr <- read.tree(text = write.tree(tr)) # otherwise, node rotation doesn't change structure


M_FIT <- readRDS("pgls-model-fits.rds")

snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']

width <- 0.1

# ----------------------------------------------------


nPreserved <- 250

preservedTips <- seq(1, Ntip(tr), by = round(Ntip(tr) / nPreserved))
preservedTips <- rev(tr$tip.label)[preservedTips]
droppedTips <- setdiff(tr$tip.label, preservedTips)

subtree <- keep.tip(tr, preservedTips)

# subtree <- ladderize(subtree, right = vertical)

subtree <- read.tree(text = write.tree(subtree)) # otherwise, node rotation doesn't change structure

preservedTips <- subtree$tip.label


tipList <- lapply(preservedTips, function(x) x)
names(tipList) <- preservedTips
patdist <- cophenetic.phylo(tr)
tol <- 1e-6
for (i in 1:length(droppedTips)) {
	xx <- patdist[droppedTips[i], preservedTips]
	assoc <- names(xx)[which(abs(xx - min(xx)) < tol)]
	# if (length(assoc) > 1) stop()
	for (j in 1:length(assoc)) {
		tipList[[assoc[j]]] <- c(tipList[[assoc[j]]], droppedTips[i])
	}
}

# sort(lengths(tipList))


# calculate weighted means

selectedTraits <- c('meanImputedCLADS', 'bamm')

wMeanList <- vector('list', length(selectedTraits))
names(wMeanList) <- selectedTraits
ptSizeList <- vector('list', length(selectedTraits))
names(ptSizeList) <- selectedTraits
for (i in 1:length(selectedTraits)) {
	
	message('\t', i, ' -- ', unlist(selectedTraits)[i])
	
	traitdat <- setNames(alldat[, unlist(selectedTraits)[i]], alldat$treename)
	traitdat <- traitdat[!is.na(traitdat)]

	identical(subtree$tip.label, names(tipList))
	
	vec <- numeric(length(tipList))	
	
	ptSizeList[[i]] <- sapply(tipList, function(x) length(intersect(names(traitdat), x)))

	if (grepl('rate|clads|bamm', unlist(selectedTraits)[i], ignore.case = TRUE)) {
		traitdat <- log(traitdat)
		
		medianTraitDat <- median(traitdat, na.rm = TRUE)
		traitDiffs <- traitdat - medianTraitDat

		for (k in 1:length(tipList)) {
			
			weights <- 1 / setNames(as.numeric(table(unlist(tipList))[tipList[[k]]]), tipList[[k]])
			
			commonsp <- intersect(names(traitdat), tipList[[k]])
			traitDiffs2 <- traitDiffs[commonsp]
			weights <- weights[commonsp]
			identical(names(traitDiffs2), names(weights))
			
			vec[k] <- sum(traitDiffs2 * weights) / sum(weights)
		}		
	} else {
	
		for (k in 1:length(tipList)) {
			
			weights <- 1 / setNames(as.numeric(table(unlist(tipList))[tipList[[k]]]), tipList[[k]])
							
			commonsp <- intersect(names(traitdat), tipList[[k]])
			traitdat2 <- traitdat[commonsp]
			weights <- weights[commonsp]
			identical(names(traitdat2), names(weights))
			# if (length(unique(traitdat2)) > 1) stop()

			
			vec[k] <- sum(traitdat2 * weights) / sum(weights)
		}
	}
	wMeanList[[i]] <- vec
}

lapply(wMeanList, head)

for (i in 1:length(ptSizeList)) {
	ptSizeList[[i]][ptSizeList[[i]] == 0] <- NA
}
ptSizeListRescaled <- lapply(ptSizeList, function(x) rescale(x, to = c(0.7, 2.2), from = range(unlist(ptSizeList), na.rm = TRUE)))


rescalePoints <- FALSE


## VERTICAL


pdf(paste0('~/Downloads/figure4_reducedTree_speciationRates.pdf'), width = 3, height = 5)

# quartz(width = 3,height=5)

mat <- matrix(c(1,1,1,2,3,4), nrow = 3, ncol = 2, byrow = FALSE)
layout(mat, widths = c(0.65, 0.35))

par(xpd = NA, mar = c(0,0,0,9), oma = c(1,1,1,1))

subtree$root.edge <- 5

qq <- plot.phylo(subtree, show.tip.label = FALSE, edge.width = 0.5, edge.color = 'black', direction = 'right', no.margin = FALSE, root.edge = TRUE)


int <- 350

refW <- 90

colpalDiv <- colorRampPalette(rev(c('#d7191c', '#fdae61','#ffffbf', '#abd9e9', 'darkblue')), alpha = TRUE)


snakeBounds <- range(which(subtree$tip.label %in% intersect(alldat[alldat$clade_Serpentes == 1, 'treename'], subtree$tip.label)))

colubriformBounds <- range(which(subtree$tip.label %in% intersect(alldat[alldat$clade_Colubriformes == 1, 'treename'], subtree$tip.label)))

sep <- grconvertX(1, from = 'user')

rect(ytop = snakeBounds[1] - 1, ybottom = snakeBounds[2] + 0.5, xleft = 230, xright = int[length(int)] + refW+15, xpd = NA, lwd = 0.5, lty = 3)



text(x = int[length(int)] + refW+5, y = 70, 'Snakes', pos = 4, xpd = NA, cex = 0.75, srt = 270)


# ClaDS
i <- 1	
			
	traitdat <- wMeanList[[i]]

	ref <- int[1]

	refVal <- median(traitdat, na.rm = TRUE)
	refVal <- rescale(refVal, from = range(traitdat, na.rm = TRUE), to = c(ref-refW, ref+refW))
	datRange <- range(traitdat, na.rm = TRUE)
	
	segments(x0 = refVal, x1 = refVal, y0 = -3, y1 = Ntip(subtree), xpd = NA, lend = 1, lwd = width, lty = 1, col = gray(0.5))
	text(x = refVal, y = -1, 'median', pos = 1, cex = 0.35)
	
	brks <- classIntervals(traitdat, n = 65, style = 'equal')

	if (grepl('rate|clads|bamm', selectedTraits[i], ignore.case = TRUE)) {
		cols <- findColours(brks, pal = colpalDiv(length(brks$brks) - 1))

		# lollipops
		segments(x0 = refVal, x1 = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), y0 = 1:Ntip(subtree), y1 = 1:Ntip(subtree), lwd = 0.15, xpd = NA, lend = 1, col = gray(0.73))
	} else {
		cols <- findColours(brks, pal = colpal(length(brks$brks) - 1))
	}
	
	if (rescalePoints) {
		ptSizes <- ptSizeListRescaled[[i]]
	} else {
		ptSizes <- 1
	}
	points(x = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), y = 1:Ntip(subtree), xpd = NA, pch = 21, bg = cols, col = 'black', lwd = 0.15, cex = ptSizes)



# dropping Iguania, Gekkonidae, Scincidae
clades <- c('Serpentes', 'lizards', 'Colubriformes', 'Hydrophis', 'Anolis', 'Liolaemus', 'Oz_Sphenomorphines', 'limb_reduced_lizards')

cladeLabels <- c('all snakes', 'all lizards', 'Colubriformes', 'Hydrophis', 'Anolis', 'Liolaemus', 'Australian\nsphenomorphines', 'limb-reduced\nnon-snakes')

fontvec <- c(1,1,1,3,3,3,1,1)

cols <- c('coral2', 'slateblue1', 'coral2', 'coral2', 'slateblue1', 'slateblue1', 'slateblue1', 'slateblue1')

par(mar = c(5,0,2,0))

plot.new()
plot.window(xlim = c(0.5, length(clades)), ylim = range(log(alldat$meanImputedCLADS), na.rm = TRUE))
# axis(1, at = 1:length(clades), labels = rep('', length(clades)), lwd = 0, lwd.ticks = 0.5, tck = -0.02)
axis(2, lwd = 0.5, lwd.ticks = 0.5, las = 1, cex.axis = 0.5, mgp = c(3, 0.5, 0), tck = -0.03)
# box(which = "plot", bty = "l", lwd = 0.5)


text(x = 1:6 + 0.6, y = -4.4, cladeLabels[1:6], srt = 50, pos = 2, xpd = NA, cex = 0.5, font = fontvec[1:6])

text(x = 7 + 0.55, y = -4.4, 'Australian', srt = 50, pos = 2, xpd = NA, cex = 0.5, font = fontvec[7])
text(x = 7 + 0.95, y = -4.4, 'sphenomorphines', srt = 50, pos = 2, xpd = NA, cex = 0.5, font = fontvec[7])

text(x = 8 + 0.55, y = -4.4, 'limb-reduced', srt = 50, pos = 2, xpd = NA, cex = 0.5, font = fontvec[8])
text(x = 8 + 0.95, y = -4.4, 'non-snakes', srt = 50, pos = 2, xpd = NA, cex = 0.5, font = fontvec[8])


mtext('log CLaDS\nspeciation rate', side = 2, line = 1.4, cex = 0.45)

for (i in 1:length(clades)) {
	
	if (paste0('clade_', clades[i]) %in% colnames(alldat)) {
		sp <- alldat[alldat[, paste0('clade_', clades[i])] == 1, 'treename']
	} else if (clades[i] == 'lizards') {
		sp <- alldat[alldat[, 'clade_Serpentes'] == 0, 'treename']
	} else if (clades[i] %in% alldat$family) {
		sp <- alldat[alldat$family == clades[i], 'treename']
	} else if (clades[i] %in% alldat$genus) {
		sp <- alldat[alldat$genus == clades[i], 'treename']
	} else if (clades[i] == 'Oz_Sphenomorphines') {
		sp <- alldat[which(as.character(alldat$geogRadiation) == 'spheno-oz'), 'treename']
	} else if (clades[i] == 'limb_reduced_lizards') {
		sp <- alldat[alldat$numberLimbs < 4 | alldat$numberDigits < 10, 'treename']
		sp <- setdiff(sp, alldat[alldat[, 'clade_Serpentes'] == 1, 'treename'])
	}
	
	sp <- intersect(sp, tr$tip.label)
	
	xx <- log(setNames(alldat[, 'meanImputedCLADS'], alldat$treename)[sp])
	
	message('\t', cladeLabels[i], ': ', round(mean(exp(xx)), 2))
	
	violin(xx, center = i, truncate = TRUE, width = 0.3, col = adjustcolor(cols[i], alpha.f = 0.75), border = NA, xpd = NA, ptCex = 0.5, barlwd = 0.5)
}


par(mar = c(2,0,2,0))


# Elongation index vs speciation rate

dat <- alldat[, c('elongationIndex', 'meanImputedCLADS')]
rownames(dat) <- alldat$treename
dat <- log(dat)
dat <- dat[complete.cases(dat), ]

plot.new()
plot.window(xlim = range(dat[,1]), ylim = range(dat[,2]))
axis(1, lwd = 0, lwd.ticks = 0.5, cex.axis = 0.5, mgp = c(3, 0.1, 0), tck = -0.02)
axis(2, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.5, mgp = c(3, 0.5, 0), tck = -0.02)
# axis(2, at = axTicks(2), labels = round(exp(axTicks(2)), 2), lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l", lwd = 0.5)
mtext('log elongation index', side = 1, line = 1, cex = 0.45)
mtext('log CLaDS\nspeciation rate', side = 2, line = 1.4, cex = 0.45)

points(dat[rownames(dat) %in% snakes,], bg = adjustcolor('coral2', alpha.f = 0.75), col = adjustcolor('black', alpha.f = 0.75), lwd = 0.45, pch = 21, cex = 0.75)
points(dat[!rownames(dat) %in% snakes,], bg = adjustcolor('slateblue1', alpha.f = 0.75), col = adjustcolor('black', alpha.f = 0.75), lwd = 0.45, pch = 21, cex = 0.75)



# ----------------------------------
# regression results
## barplot of R2 values

par(mar = c(5,0,2,0))

r2vals <- sapply(M_FIT, function(x) summary(x)$adj.r.squared)
r2vals[r2vals < 0] <- 0 # convert negative adjusted R2 to zero
sort(round(r2vals, 3))

# better names for figure
newnames <- c(M_TROPHIC1 = 'diet', M_TROPHIC2 = 'foraging mode', M_TROPHIC3 = 'chemosensory', M_PARITY = 'parity mode', M_KINESIS = 'skull kinesis', M_BIOCLIM = 'climate', M_MORPH = 'body size/shape', M_MERISTIC1 = 'digit/limb count', M_MERISTIC2 = 'vertebral count', M_SKULL = 'skull shape', M_GEO = 'geography', M_BIOGEO = 'biogeog. theater', M_CLADE = 'clade membership')

names(r2vals) <- newnames[names(r2vals)]

M_FIT2 <- M_FIT[order(r2vals, decreasing = TRUE)]
names(M_FIT2) <- newnames[names(M_FIT2)]

r2vals <- sort(r2vals, decreasing = TRUE)


xx <- barplot(r2vals, ylim = c(0,1), las = 3, names.arg = NA, axes = FALSE, lwd = 0.5, lend = 1)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), lwd = 0.5, lwd.ticks = 0.5, las = 1, cex.axis = 0.5, mgp = c(3, 0.5, 0), tck = -0.02)
text(x = xx[,1] + 1, y = -0.04, names(r2vals), srt = 50, pos = 2, xpd = NA, cex = 0.5)
mtext(expression(paste('phylogenetic R' ^2)), side = 2, line = 1.4, cex = 0.45)


text(x = grconvertX(0.05, from = 'ndc', to = 'user'), y = grconvertY(0.96, from = 'ndc', to = 'user'), 'A.', font = 2, cex = 1)
text(x = grconvertX(0.55, from = 'ndc', to = 'user'), y = grconvertY(0.96, from = 'ndc', to = 'user'), 'B.', font = 2, cex = 1)
text(x = grconvertX(0.55, from = 'ndc', to = 'user'), y = grconvertY(0.65, from = 'ndc', to = 'user'), 'C.', font = 2, cex = 1)
text(x = grconvertX(0.55, from = 'ndc', to = 'user'), y = grconvertY(0.33, from = 'ndc', to = 'user'), 'D.', font = 2, cex = 1)



dev.off()

