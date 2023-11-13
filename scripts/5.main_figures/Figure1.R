# Figure 1:


setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

require(ape)
require(classInt)
require(phytools)
require(scales)

alldat <- read.csv('./data/alldat.csv')

# nodeRotationOption <- 1 # good for vertical
nodeRotationOption <- 2 # good for horizontal

tr <- read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
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



snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']

width <- 0.1




# ------------------------------------------

selectedTraits <- 	list(
						"Rates" = list(
							climate = 'climateNicheRate',
							'vertebral' = 'vertebralRate',
							"skull" = 'skullRate',
							"elong." = 'elongationRate',
							'diet' = 'dietRate'
						),
						"Net innovation" = list(
							'skull' = 'skullAncDist',
							'vert' = 'ancDistVert',
							'elong' = 'ancDistElongationIndex',
							'diet' = 'ancDistDiet'
						)
					)

table(unlist(selectedTraits) %in% colnames(alldat))



# ----------------------------------------------


# for nPreserved species, get pruned tree with just those taxa
# for taxa that are left out, associate them with the most closely related preserved tips.
# if multiple preserved tips are equally closely related, associate with each of them.

# convert chemosensory index to factor
# alldat[, 'chemosensory_index'] <- as.factor(alldat[, 'chemosensory_index'])

nPreserved <- 250

preservedTips <- seq(1, Ntip(tr), by = round(Ntip(tr) / nPreserved))
preservedTips <- tr$tip.label[preservedTips]
droppedTips <- setdiff(tr$tip.label, preservedTips)

subtree <- keep.tip(tr, preservedTips)

# subtree <- ladderize(subtree, right = FALSE)


subtree <- read.tree(text = write.tree(subtree)) # otherwise, node rotation doesn't change structure

identical(preservedTips, subtree$tip.label)

preservedTips <- subtree$tip.label

plot.phylo(subtree, show.tip.label = FALSE, direction = 'upwards')
segments(y0 = 215, y1 = 215, x0 = min(which(subtree$tip.label %in% snakes)), x1 = max(which(subtree$tip.label %in% snakes)), lend = 1, col = 'red', xpd = NA, lwd = 10)
points(y = rep(215, length(grep('Liolaemus', subtree$tip.label, value = TRUE))), x = grep('Liolaemus', subtree$tip.label), pch = 23, xpd = NA)
points(y = rep(215, length(grep('Bradypodion|Trioceros', subtree$tip.label, value = TRUE))), x = grep('Bradypodion|Trioceros', subtree$tip.label), pch = 25, xpd = NA)


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
wMeanList <- vector('list', length(unlist(selectedTraits)))
names(wMeanList) <- unlist(selectedTraits)
wMeanDiffList <- vector('list', length(unlist(selectedTraits)))
names(wMeanDiffList) <- unlist(selectedTraits)
ptSizeList <- vector('list', length(unlist(selectedTraits)))
names(ptSizeList) <- unlist(selectedTraits)
for (i in 1:length(unlist(selectedTraits))) {
	
	message('\t', i, ' -- ', unlist(selectedTraits)[i])
	
	traitdat <- setNames(alldat[, unlist(selectedTraits)[i]], alldat$treename)
	traitdat <- traitdat[!is.na(traitdat)]

	identical(subtree$tip.label, names(tipList))
	
	vec <- numeric(length(tipList))
	vecDiff <- numeric(length(tipList))
	
	ptSizeList[[i]] <- sapply(tipList, function(x) length(intersect(names(traitdat), x)))
	
	# for categorical data, just preserve majority value
	if (is.factor(traitdat)) {
		traitdat <- setNames(as.character(traitdat), names(traitdat))
		
		# for parity, for each tip, calc percent of taxa that are viviparous
		# for chemosensory, calc average chemosensory index <- therefore it gets treated as a continuous variable
		
		if (unlist(selectedTraits)[i] == 'parity') {
			
			for (k in 1:length(tipList)) {
				
				commonsp <- intersect(names(traitdat), tipList[[k]])
				vec[k] <- sum(traitdat[commonsp] != 'Oviparous') / length(commonsp)
				
			}
		}
				
	} else {			

		if (grepl('rate|clads', unlist(selectedTraits)[i], ignore.case = TRUE)) {
			traitdat <- log(traitdat)
		}
		
		medianTraitDat <- median(traitdat, na.rm = TRUE)
		traitDiffs <- traitdat - medianTraitDat
	
		for (k in 1:length(tipList)) {
			
			weights <- 1 / setNames(as.numeric(table(unlist(tipList))[tipList[[k]]]), tipList[[k]])
			
			commonsp <- intersect(names(traitdat), tipList[[k]])
			traitDiffs2 <- traitDiffs[commonsp]
			traitdat2 <- traitdat[commonsp]
			weights <- weights[commonsp]
			identical(names(traitDiffs2), names(weights))
			
			vec[k] <- sum(traitdat2 * weights) / sum(weights)
			vecDiff[k] <- sum(traitDiffs2 * weights) / sum(weights)
		}
	}		

	wMeanList[[i]] <- vec
	wMeanDiffList[[i]] <- vecDiff
}

lapply(wMeanList, head)

for (i in 1:length(ptSizeList)) {
	ptSizeList[[i]][ptSizeList[[i]] == 0] <- NA
}
ptSizeListRescaled <- lapply(ptSizeList, function(x) rescale(x, to = c(0.7, 2.2), from = range(unlist(ptSizeList), na.rm = TRUE)))


## 
ResizePoints <- FALSE

style <- 'points'

pdf(paste0('~/Downloads/figure1_reducedTree_', style, '.pdf'), width = 12, height = 12)

mat <- matrix(c(4,4,1,1,2,3), nrow = 3, ncol = 2, byrow = TRUE)
mat <- matrix(c(2,3,4,4,1,1), nrow = 3, ncol = 2, byrow = TRUE)
layout(mat, heights = c(0.25, 0.60, 0.15))

par(xpd = NA, oma = c(1,3,1,1))

subtree$root.edge <- 5

qq <- plot.phylo(subtree, show.tip.label = FALSE, edge.width = 0.5, edge.color = 'black', direction = 'upwards', no.margin = TRUE, root.edge = TRUE)

# axisPhylo(side = 2, line = -2.5, cex.axis = 0.75)
oldest <- max(branching.times(subtree)) + 5
segments(x0 = -2, x1 = -2, y0 = 0, y1 = oldest, xpd = NA, lend = 1, lwd = 0.5)
tx <- seq(from = 0, to = 200, by = 50)
segments(x0 = -3, x1 = -2, y0 = oldest - tx, y1 = oldest - tx, xpd = NA, lend = 1, lwd = 0.5)
text(x = -2.5, y = oldest - tx, labels = tx, pos = 2, cex = 0.5)

# label major nodes
toLabel <- c('Serpentes', 'Alethinophidia', 'Colubriformes', 'Gekkota', 'Scincoidea', 'Anguimorpha', 'Iguania', 'Amphisbaenia', 'Lacertoidea')
adjList <- list(
	c(-0.1, 1.1),
	c(-0.1, 1.1),
	c(-0.1, 1.1),
	c(-0.1, 1.1),
	c(-0.1, 1.1),
	c(-0.1, 1.1),
	c(0, -1), # iguania
	c(-0.1, 1.1), # amphisbaenia
	c(1.1, 1.1)  # Lacertoidea
)
cl <- grep('clade_', colnames(alldat), value = TRUE)
for (i in 1:length(toLabel)) {
	nn <- getMRCA(subtree, intersect(subtree$tip.label, alldat[alldat[, paste0('clade_', toLabel[i])] == 1, 'treename']))
	if (length(nn) > 0) {
		nodelabels(toLabel[i], node = nn, frame = 'none', cex = 1, adj = adjList[[i]])
		nodelabels(node = nn, pch = 21, cex = 0.5, bg = 'black')
	}
}

# ----------------------
maxPlot <- 1100
int <- seq(260, 950, length.out = length(unlist(selectedTraits)))


refW <- diff(int)[1]/2 * 0.85
if ((int[length(int)] + refW) > maxPlot) {
	refW <- maxPlot - int[length(int)]
}

for (i in 2:length(selectedTraits)) {
	ii <- (which(unlist(selectedTraits) == selectedTraits[[i]][1])) : length(int)
	int[ii] <- int[ii] + 30
	segments(x0 = grconvertX(0, from = 'nfc', to ='user'), x1 = Ntip(subtree), y0 = mean(c(int[ii[1]], int[ii[1] - 1])), y1 = mean(c(int[ii[1]], int[ii[1] - 1])), lwd = 0.4)
}

width2 <- 0.4

# axis(2, at = seq(0, 1200, by = 50), xpd = NA, line = -1)
# rect(xleft = 1, ybottom = int - refW, xright = Ntip(subtree), ytop = int + refW, xpd = NA)

rect(xleft = 1, ybottom = int - refW, xright = Ntip(subtree), ytop = int + refW, xpd = NA, col = gray(0.95), border = NA)


colpal <- viridis::plasma

colpalDiv <- colorRampPalette(rev(c('#d7191c', '#fdae61','#ffffbf', '#abd9e9', 'darkblue')), alpha = TRUE)

segments(
	x0 = rescale(1:100, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	x1 = rescale(1:100, , from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y0 = rescale(rep(0, 99), from = c(0,1), to = c(100, 115)), 
	y1 = rescale(rep(1, 99), from = c(0,1), to = c(100, 115)), 
	col = colpalDiv(100), lend = 1, lwd = 2)
segments(	
	x0 = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	x1 = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y0 = rescale(0, from = c(0,1), to = c(100, 115)), 
	y1 = rescale(1.4, from = c(0,1), to = c(100, 115)), 
	col = 'black', lend = 1, lwd = 0.5)
text(
	x = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(1.4, from = c(0,1), to = c(100, 115)),
	labels = 'median', pos = 3, cex = 0.75)
text(
	x = rescale(1, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(0.5, from = c(0,1), to = c(100, 115)),
	labels = 'low', pos = 2, cex = 1)
text(
	x = rescale(100, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(0.5, from = c(0,1), to = c(100, 115)),
	labels = 'high', pos = 4, cex = 1)
	
	
# scaling legend
points(x = grconvertX(seq(0.85, 0.92, length.out = 4), from='npc', to = 'user'), y = rep(70, 4), cex = seq(min(unlist(ptSizeListRescaled), na.rm = TRUE), max(unlist(ptSizeListRescaled), na.rm = TRUE), length.out = 4))
text(x=grconvertX(seq(0.85, 0.92, length.out = 4), from='npc', to = 'user'), y = rep(67, 4), labels = seq(min(unlist(ptSizeList), na.rm = TRUE), max(unlist(ptSizeList), na.rm = TRUE), length.out = 4), pos = 1)
text(x=grconvertX(0.945, from='npc', to = 'user'), y = 67, labels = 'sp.', pos = 1)

snakeBounds <- range(which(subtree$tip.label %in% intersect(alldat[alldat$clade_Serpentes == 1, 'treename'], subtree$tip.label)))

colubriformBounds <- range(which(subtree$tip.label %in% intersect(alldat[alldat$clade_Colubriformes == 1, 'treename'], subtree$tip.label)))

sep <- grconvertX(1, from = 'user')

rect(xleft = snakeBounds[1] - 0.5, xright = snakeBounds[2] + 1, ybottom = 220, ytop = int[length(int)] + refW+30, xpd = NA, lwd = 0.5, lty = 3)

rect(xleft = colubriformBounds[1] - 0.5, xright = colubriformBounds[2] + 1, ybottom = 220, ytop = int[length(int)] + refW + 12, xpd = NA, lwd = 0.5, lty = 3)


text(x = grconvertX(0.8, from = 'nfc', to ='user'), y = int[length(int)] + refW + 19, 'Colubriformes', pos = 4, xpd = NA, cex = 1.2)
text(x = grconvertX(0.9, from = 'nfc', to ='user'), y = int[length(int)] + refW + 37, 'Serpentes', pos = 4, xpd = NA, cex = 1.2)

refW <- refW * 0.90

for (i in 1:length(selectedTraits)) {
	
	for (j in 1:length(selectedTraits[[i]])) {
			
		ij <- which(names(wMeanList) == selectedTraits[[i]][[j]])
		ref <- int[ij]
		
			traitdat <- wMeanDiffList[[ij]]

			refVal <- median(traitdat, na.rm = TRUE)
			refVal <- rescale(refVal, from = range(traitdat, na.rm = TRUE), to = c(ref-refW, ref+refW))
			datRange <- range(traitdat, na.rm = TRUE)

			segments(x0 = 1, x1 = Ntip(subtree), y0 = refVal, y1 = refVal, xpd = NA, lend = 1, lwd = width, lty = 1, col = gray(0.5))
			
			brks <- classIntervals(traitdat, n = 65, style = 'equal')
			cols <- findColours(brks, pal = colpalDiv(length(brks$brks) - 1))

			# lollipops
			segments(x0 = 1:Ntip(subtree), x1 = 1:Ntip(subtree), y0 = refVal, y1 = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), lwd = 0.15, xpd = NA, lend = 1)
					
		if (ResizePoints) {
			ptSizes <- ptSizeListRescaled[[ij]]
		} else {
			ptSizes <- 1.3
		}
		
		points(x = 1:Ntip(subtree), y = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), xpd = NA, pch = 21, bg = cols, col = 'black', lwd = 0.15, cex = ptSizes)

			
		if (grepl('rate|clads', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			lab <- bquote(TR[.(names(selectedTraits[[i]])[j])])
		} else if (grepl('ancDist|chemo', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			# lab <- bquote(Phi[.(names(selectedTraits[[i]])[j])])
			lab <- bquote(psi[.(names(selectedTraits[[i]])[j])])
		} else {
			lab <- names(selectedTraits[[i]])[j]
		}
		text(x = grconvertX(0.03, from = 'nfc', to ='user'), y = ref, labels = lab, xpd = NA, pos = 2, srt = 0, cex = 1)
		
		if (j == 1) {
			text(x = grconvertX(-0.01, from = 'nfc', to ='user'), y = mean(int[unlist(selectedTraits) %in% selectedTraits[[i]]]), labels = names(selectedTraits)[i], xpd = NA, pos = 3, srt = 90, cex = 1.5, font = 2)
		}
	}	
}


dev.off()



