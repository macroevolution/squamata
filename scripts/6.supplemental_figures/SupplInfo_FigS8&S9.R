setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

require(ape)
require(classInt)
require(phytools)
require(scales)

alldat <- read.csv('./data/alldat.csv')

# nodeRotationOption <- 1 # good for vertical
nodeRotationOption <- 2 # good for horizontal

tr <- read.tree("././data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
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




nPreserved <- 250

preservedTips <- seq(1, Ntip(tr), by = round(Ntip(tr) / nPreserved))
preservedTips <- tr$tip.label[preservedTips]
droppedTips <- setdiff(tr$tip.label, preservedTips)

subtree <- keep.tip(tr, preservedTips)



subtree <- read.tree(text = write.tree(subtree)) # otherwise, node rotation doesn't change structure

identical(preservedTips, subtree$tip.label)

preservedTips <- subtree$tip.label

# plot.phylo(subtree, show.tip.label = FALSE, direction = 'upwards')
# segments(y0 = 215, y1 = 215, x0 = min(which(subtree$tip.label %in% snakes)), x1 = max(which(subtree$tip.label %in% snakes)), lend = 1, col = 'red', xpd = NA, lwd = 10)
# points(y = rep(215, length(grep('Liolaemus', subtree$tip.label, value = TRUE))), x = grep('Liolaemus', subtree$tip.label), pch = 23, xpd = NA)
# points(y = rep(215, length(grep('Bradypodion|Trioceros', subtree$tip.label, value = TRUE))), x = grep('Bradypodion|Trioceros', subtree$tip.label), pch = 25, xpd = NA)


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

# convert kinesis to ordered factor to ensure proper plotting
alldat$combinedKinesis <- factor(alldat$combinedKinesis, levels = c("akinetic", "midKinesis", "lowKinesis", "highKinesis", "hyperkinesis"), ordered = TRUE)

# and convert to numeric for the purposes of plotting
alldat$combinedKinesis <- as.numeric(alldat$combinedKinesis)


###########################################################

## VERSION ORGANIZED AROUND PRIMARY/SECONDARY/TYPE 1/2 TRAITS


selectedTraits1 <- list(
			'morphological' = list(
				'log mass' = 'mass',
				'log SVL' = 'completeSVL',
				'elongation\nindex' = 'elongationIndex',
				'# presacral\nvertebrae' = 'numberPresacralVert',
				'vertInnovation' = 'ancDistVert',
				'vertRate' = 'vertebralRate',
				'skull\nshape\n(PC1)' = 'skullPC1'
			),
			'environmental' = list(
				'climate\n(PC1)' = 'climPC1',
				'climateRate' = 'climateNicheRate'
			),
			'ecological' = list(
				'chemosensory\ninnovation' = 'chemosensory_index',
				'diet (PC1)' = 'dietPC1',
				'diet\nbreadth' = 'dietBreadth'
			)
)

selectedTraits2 <- list(
			'morphological' = list(
				'# digits' = 'numberDigits',
				'# limbs' = 'numberLimbs',
				'skull\nkinesis' = 'combinedKinesis',	
				'prehension\nmechanism' = 'prehensionMechanism'
			),
			'ecological' = list(	
				'parity' = 'parity',
				'foraging\nmode' = 'foragingMode'
			),
			'Biogeographic traits' = list(
				'latitude' = 'centroidLat',
				'elevation' = 'elev'
			)
)

# we reverse it because the plotting plots from the bottom up
selectedTraits1 <- lapply(selectedTraits1, rev)
selectedTraits1 <- rev(selectedTraits1)

selectedTraits2 <- lapply(selectedTraits2, rev)
selectedTraits2 <- rev(selectedTraits2)


discrete <- c('parity', 'prehensionMechanism', 'foragingMode')




# calculate weighted means
wMeanList <- vector('list', length(unlist(selectedTraits1)))
names(wMeanList) <- unlist(selectedTraits1)
wMeanDiffList <- vector('list', length(unlist(selectedTraits1)))
names(wMeanDiffList) <- unlist(selectedTraits1)
ptSizeList <- vector('list', length(unlist(selectedTraits1)))
names(ptSizeList) <- unlist(selectedTraits1)
for (i in 1:length(unlist(selectedTraits1))) {
	
	message('\t', i, ' -- ', unlist(selectedTraits1)[i])


	if (unlist(selectedTraits1)[i] == 'kinesis') {
		
		traitdat <- rep(NA, Ntip(tr))
		names(traitdat) <- tr$tip.label
		commonsp <- intersect(names(kinesis2), tr$tip.label)
		traitdat <- traitdat[commonsp]
		traitdat[match(names(kinesis2[commonsp]), names(traitdat))] <- kinesis2[commonsp]

	} else {

		traitdat <- setNames(alldat[, unlist(selectedTraits1)[i]], alldat$treename)
		# range(traitdat, na.rm = TRUE)

	}
		
	traitdat <- traitdat[!is.na(traitdat)]

	identical(subtree$tip.label, names(tipList))
	
	vec <- numeric(length(tipList))
	vecDiff <- numeric(length(tipList))
	
	ptSizeList[[i]] <- sapply(tipList, function(x) length(intersect(names(traitdat), x)))
	
	# for categorical data, just preserve majority value
	if (unlist(selectedTraits1)[i] %in% discrete) {
		traitdat <- setNames(as.character(traitdat), names(traitdat))
		
		# for parity, for each tip, calc percent of taxa that are viviparous
		# for chemosensory, calc average chemosensory index <- therefore it gets treated as a continuous variable
		
		if (unlist(selectedTraits1)[i] == 'parity') {
			
			for (k in 1:length(tipList)) {
				
				commonsp <- intersect(names(traitdat), tipList[[k]])
				vec[k] <- sum(traitdat[commonsp] != 'Oviparous') / length(commonsp)
				
			}
		}

		if (unlist(selectedTraits1)[i] == 'foragingMode') {
			
			for (k in 1:length(tipList)) {
				
				commonsp <- intersect(names(traitdat), tipList[[k]])
				vec[k] <- sum(traitdat[commonsp] != 'ambush') / length(commonsp)
				
			}
		}

		if (unlist(selectedTraits1)[i] == 'kinesis') {
			
			for (k in 1:length(tipList)) {
				
				commonsp <- intersect(names(traitdat), tipList[[k]])
				vec[k] <- sum(traitdat[commonsp] != 'akinetic') / length(commonsp)
				
			}
		}

		if (unlist(selectedTraits1)[i] == 'prehensionMechanism') {
			
			for (k in 1:length(tipList)) {
				
				commonsp <- intersect(names(traitdat), tipList[[k]])
				vec[k] <- sum(traitdat[commonsp] %in% c('lingual', 'both')) / length(commonsp)
				
			}
		}
				
	} else {

		if (grepl('rate|clads|svl|mass', unlist(selectedTraits1)[i], ignore.case = TRUE)) {
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
ptSizeListRescaled <- lapply(ptSizeList, function(x) rescale(x, to = c(0.5, 2), from = range(unlist(ptSizeList), na.rm = TRUE)))

style <- 'points'

ResizePoints <- FALSE


## FIRST PART WITH PRIMARY TRAITS ONLY


pdf(paste0('~/Downloads/figS8.pdf'), width = 11, height = 10)

mat <- matrix(c(2,1), nrow = 2, ncol = 1, byrow = TRUE)
layout(mat, heights = c(0.85, 0.15))

par(xpd = NA, oma = c(2,8,3,1))

subtree$root.edge <- 5

qq <- plot.phylo(subtree, show.tip.label = FALSE, edge.width = 0.5, edge.color = 'black', direction = 'upwards', no.margin = TRUE, root.edge = TRUE)

# axisPhylo(side = 2, line = -1, cex.axis = 0.75)
oldest <- max(branching.times(subtree)) + 5
segments(x0 = -2, x1 = -2, y0 = 0, y1 = oldest, xpd = NA, lend = 1, lwd = 0.5)
tx <- seq(from = 0, to = 200, by = 50)
segments(x0 = -3, x1 = -2, y0 = oldest - tx, y1 = oldest - tx, xpd = NA, lend = 1, lwd = 0.5)
text(x = -2.5, y = oldest - tx, labels = tx, pos = 2, cex = 0.5)


# # 
# # label major nodes
# toLabel <- c('Serpentes', 'Alethinophidia', 'Colubriformes', 'Gekkota', 'Scincoidea', 'Teioidea', 'Anguiformes', 'Iguania', 'Amphisbaenia', 'Lacertoidea')
# adjList <- list(
	# c(-0.1, 1.1),
	# c(-0.1, 1.1),
	# c(-0.1, 1.1),
	# c(-0.1, 1.1),
	# c(1.1, 1.1),
	# c(1.1, 1.1),
	# c(-0.5, 1), #Anguiformes
	# c(0, -0.9), # iguania
	# c(-0.05, 1.1), # amphisbaenia
	# c(1.3, 0.9)  # Lacertoidea
# )
# cl <- grep('clade_', colnames(alldat), value = TRUE)
# for (i in 1:length(toLabel)) {
	# nn <- getMRCA(subtree, intersect(subtree$tip.label, alldat[alldat[, paste0('clade_', toLabel[i])] == 1, 'treename']))
	# if (length(nn) > 0) {
		# nodelabels(toLabel[i], node = nn, frame = 'none', cex = 0.9, adj = adjList[[i]])
		# nodelabels(node = nn, pch = 21, cex = 0.5, bg = 'black')
	# }
# }


nodeLabels <- c(A = 'Dibamidae', B = 'Gekkota', C = 'Scincoidea', D = 'Teioidea', E = 'Amphisbaenia', F = 'Lacertoidea', G = 'Anguiformes', H = 'Iguania', I = 'Snakes', J = 'Scolecophidia', K = 'Alethinophidia', L = 'Colubriformes')
points(x = 1, y = 50, pch = 21, bg = 'white', cex = 1, lwd = 0.5)
text(x = 1, y = 50, 'A', cex = 0.3)
for (i in 2:length(nodeLabels)) {
	cc <- paste0('clade_', nodeLabels[i])
	if (nodeLabels[i] == 'Snakes') cc <- 'clade_Serpentes'
	nn <- getMRCA(subtree, intersect(subtree$tip.label, alldat[alldat[, cc] == 1, 'treename']))
	if (length(nn) > 0) {
		nodelabels(node = nn, pch = 21, cex = 1, bg = 'white', lwd = 0.5)
		nodelabels(names(nodeLabels)[i], node = nn, frame = 'none', cex = 0.3)
		
	}
}



# ----------------------
maxPlot <- 1500
int <- seq(300, 1400, length.out = length(unlist(selectedTraits1)))


refW <- diff(int)[1]/2 * 0.85
if ((int[length(int)] + refW) > maxPlot) {
	refW <- maxPlot - int[length(int)]
}

for (i in 2:length(selectedTraits1)) {
	ii <- (which(unlist(selectedTraits1) == selectedTraits1[[i]][1])) : length(int)
	int[ii] <- int[ii] + 30
	segments(x0 = grconvertX(0, from = 'nfc', to ='user'), x1 = Ntip(subtree), y0 = mean(c(int[ii[1]], int[ii[1] - 1])), y1 = mean(c(int[ii[1]], int[ii[1] - 1])), lwd = 0.4)
}

width2 <- 0.4

# axis(2, at = seq(0, 2000, by = 50), xpd = NA, line = -1)
# rect(xleft = 1, ybottom = int - refW, xright = Ntip(subtree), ytop = int + refW, xpd = NA)

rect(xleft = 1, ybottom = int - refW, xright = Ntip(subtree), ytop = int + refW, xpd = NA, col = gray(0.95), border = NA)


colpal <- viridis::plasma
#colpal <- viridis::inferno
# colpal <- colorRampPalette(BAMMtools:::palettes$RdYlBu)
#colpal <- viridis::plasma(5)
#colpal[1] <- 'darkblue'
#colpal <- colorRampPalette(colpal)

colpalDiv <- colorRampPalette(rev(c('#d7191c', '#fdae61','#ffffbf', '#abd9e9', 'darkblue')), alpha = TRUE)

# plot(1:20, rep(1, 20), bg = colpal(20), pch = 21, cex = 5)
# points(1:20, rep(0.8, 20), bg = colpalDiv(20), pch = 21, cex = 5)

segments(
	x0 = rescale(1:100, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	x1 = rescale(1:100, , from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y0 = rescale(rep(0, 99), from = c(0,1), to = c(50, 65)), 
	y1 = rescale(rep(1, 99), from = c(0,1), to = c(50, 65)), 
	col = colpalDiv(100), lend = 1, lwd = 2)
segments(	
	x0 = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	x1 = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y0 = rescale(0, from = c(0,1), to = c(50, 65)), 
	y1 = rescale(1.4, from = c(0,1), to = c(50, 65)), 
	col = 'black', lend = 1, lwd = 0.5)
text(
	x = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(1, from = c(0,1), to = c(50, 65)),
	labels = 'median', pos = 3, cex = 0.75)
text(
	x = rescale(1, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(0.5, from = c(0,1), to = c(50, 65)),
	labels = 'low', pos = 2, cex = 1)
text(
	x = rescale(100, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(0.5, from = c(0,1), to = c(50, 65)),
	labels = 'high', pos = 4, cex = 1)
	
# # scaling legend
# points(x = grconvertX(seq(0.85, 0.92, length.out = 4), from='npc', to = 'user'), y = rep(30, 4), cex = seq(min(unlist(ptSizeListRescaled), na.rm = TRUE), max(unlist(ptSizeListRescaled), na.rm = TRUE), length.out = 4))
# text(x=grconvertX(seq(0.85, 0.92, length.out = 4), from='npc', to = 'user'), y = rep(25, 4), labels = seq(min(unlist(ptSizeList), na.rm = TRUE), max(unlist(ptSizeList), na.rm = TRUE), length.out = 4), pos = 1, cex = 0.8)
# text(x=grconvertX(0.945, from='npc', to = 'user'), y = 25, labels = 'sp.', pos = 1, cex = 0.8)


snakeBounds <- range(which(subtree$tip.label %in% intersect(alldat[alldat$clade_Serpentes == 1, 'treename'], subtree$tip.label)))

colubriformBounds <- range(which(subtree$tip.label %in% intersect(alldat[alldat$clade_Colubriformes == 1, 'treename'], subtree$tip.label)))

sep <- grconvertX(1, from = 'user')

# rect(xleft = snakeBounds[1] - 15, xright = snakeBounds[2] + 0.5, ybottom = 240 - refW, ytop = int[length(int)] + refW+25, xpd = NA, lwd = 0.5, col = adjustcolor(gray(0.9), alpha.f = 0.5), border = NA)

# rect(xleft = colubriformBounds[1] - 10, xright = colubriformBounds[2] + 0.5, ybottom = 240 - refW, ytop = int[length(int)] + refW + 12, xpd = NA, lwd = 0.5, col = adjustcolor(gray(0.9), alpha.f = 0.5), border = NA)

rect(xleft = snakeBounds[1] - 0.5, xright = snakeBounds[2] + 1, ybottom = 240 - 15, ytop = int[length(int)] + refW+35, xpd = NA, lwd = 0.5, lty = 3)

rect(xleft = colubriformBounds[1] - 0.5, xright = colubriformBounds[2] + 1, ybottom = 240 - 15, ytop = int[length(int)] + refW + 12, xpd = NA, lwd = 0.5, lty = 3)


text(x = grconvertX(0.75, from = 'nfc', to ='user'), y = int[length(int)] + refW + 19, 'Colubriformes', pos = 4, xpd = NA, cex = 0.8)
text(x = grconvertX(0.85, from = 'nfc', to ='user'), y = int[length(int)] + refW + 42, 'Serpentes', pos = 4, xpd = NA, cex = 0.8)


refW <- refW * 0.75

#counter <- 1
for (i in 1:length(selectedTraits1)) {

	for (j in 1:length(selectedTraits1[[i]])) {
			
		ij <- which(names(wMeanList) == selectedTraits1[[i]][[j]])
		ref <- int[ij]


		if (grepl('rate|clads', selectedTraits1[[i]][[j]], ignore.case = TRUE)) {

			traitdat <- wMeanDiffList[[ij]]

			refVal <- median(traitdat, na.rm = TRUE)
			refVal <- rescale(refVal, from = range(traitdat, na.rm = TRUE), to = c(ref-refW, ref+refW))
			datRange <- range(traitdat, na.rm = TRUE)

			segments(x0 = 1, x1 = Ntip(subtree), y0 = refVal, y1 = refVal, xpd = NA, lend = 1, lwd = width, lty = 1, col = gray(0.5))
			
			brks <- classIntervals(traitdat, n = 65, style = 'equal')
			cols <- findColours(brks, pal = colpalDiv(length(brks$brks) - 1))

			# lollipops
#			segments(x0 = 1:Ntip(subtree), x1 = 1:Ntip(subtree), y0 = refVal, y1 = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), lwd = 0.15, xpd = NA, lend = 1)
			
		} else {
			
			if (selectedTraits1[[i]][[j]] %in% discrete) {

				traitdat <- wMeanList[[ij]]
				medianDiff <- traitdat
			
			} else {
				
				traitdat <- wMeanList[[ij]]
				medianDiff <- wMeanDiffList[[ij]]

			}
			
#			refVal <- ref-refW
#			datRange <- c(0, max(traitdat, na.rm = TRUE))

			refVal <- median(traitdat, na.rm = TRUE)
			refVal <- rescale(refVal, from = range(traitdat, na.rm = TRUE), to = c(ref-refW, ref+refW))
			datRange <- range(traitdat, na.rm = TRUE)

			segments(x0 = 1, x1 = Ntip(subtree), y0 = refVal, y1 = refVal, xpd = NA, lend = 1, lwd = width, lty = 1, col = gray(0.5))

			# brks <- classIntervals(medianDiff, n = 65, style = 'equal')
			brks <- classIntervals(round(medianDiff, 10), n = 65, style = 'equal')
			cols <- findColours(brks, pal = colpalDiv(length(brks$brks) - 1))
			
		}
		
		if (ResizePoints) {
			ptSizes <- ptSizeListRescaled[[ij]]
		} else {
			ptSizes <- 0.8
		}
		
#		# lollipops
#		segments(x0 = 1:Ntip(subtree), x1 = 1:Ntip(subtree), y0 = refVal, y1 = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), lwd = 0.15, xpd = NA, lend = 1)
	
		points(x = 1:Ntip(subtree), y = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), xpd = NA, pch = 21, bg = cols, col = 'black', lwd = 0.15, cex = ptSizes)

		lab <- names(selectedTraits1[[i]])[j]
		if (grepl('chemosensory', lab)) {
			lab <- bquote(psi[chem])
		
		} else if (grepl('vertInnovation', lab)) {
			lab <- 	bquote(psi[vert])

		} else if (grepl('climateRate', lab)) {
			lab <- 	bquote(TR[climate])

		} else if (grepl('vertRate', lab)) {
			lab <- 	bquote(TR[vert])
		}
		
		# if (grepl('rate|clads', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			# lab <- bquote(beta[.(names(selectedTraits[[i]])[j])])
		# } else if (grepl('ancDist|chemo', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			# lab <- bquote(Phi[.(names(selectedTraits[[i]])[j])])
		# } else {
			# lab <- names(selectedTraits[[i]])[j]
		# }
		text(x = grconvertX(0.03, from = 'nfc', to ='user'), y = ref, labels = lab, xpd = NA, pos = 2, srt = 0, cex = 0.75)
		
		if (j == 1) {
			text(x = grconvertX(-0.06, from = 'nfc', to ='user'), y = mean(int[unlist(selectedTraits1) %in% selectedTraits1[[i]]]), labels = names(selectedTraits1)[i], xpd = NA, pos = 3, srt = 90)
		}
	}	
}

text(x = grconvertX(-0.12, from = 'nfc', to ='user'), y = grconvertY(0.55, from = 'ndc'), labels = 'primary traits', xpd = NA, pos = 3, srt = 90, font = 2)

dev.off()




############################################
## SECONDARY TRAITS - FigS8



# calculate weighted means
wMeanList <- vector('list', length(unlist(selectedTraits2)))
names(wMeanList) <- unlist(selectedTraits2)
wMeanDiffList <- vector('list', length(unlist(selectedTraits2)))
names(wMeanDiffList) <- unlist(selectedTraits2)
ptSizeList <- vector('list', length(unlist(selectedTraits2)))
names(ptSizeList) <- unlist(selectedTraits2)
for (i in 1:length(unlist(selectedTraits2))) {
	
	message('\t', i, ' -- ', unlist(selectedTraits2)[i])


	# if (unlist(selectedTraits2)[i] == 'kinesis') {
		
		# traitdat <- rep(NA, Ntip(tr))
		# names(traitdat) <- tr$tip.label
		# commonsp <- intersect(names(kinesis2), tr$tip.label)
		# traitdat <- traitdat[commonsp]
		# traitdat[match(names(kinesis2[commonsp]), names(traitdat))] <- kinesis2[commonsp]

	# } else {

		traitdat <- setNames(alldat[, unlist(selectedTraits2)[i]], alldat$treename)
		# range(traitdat, na.rm = TRUE)

	# }
		
	traitdat <- traitdat[!is.na(traitdat)]

	identical(subtree$tip.label, names(tipList))
	
	vec <- numeric(length(tipList))
	vecDiff <- numeric(length(tipList))
	
	ptSizeList[[i]] <- sapply(tipList, function(x) length(intersect(names(traitdat), x)))
	
	# for categorical data, just preserve majority value
	if (unlist(selectedTraits2)[i] %in% discrete) {
		traitdat <- setNames(as.character(traitdat), names(traitdat))
		
		# for parity, for each tip, calc percent of taxa that are viviparous
		# for chemosensory, calc average chemosensory index <- therefore it gets treated as a continuous variable
		
		if (unlist(selectedTraits2)[i] == 'parity') {
			
			for (k in 1:length(tipList)) {
				
				commonsp <- intersect(names(traitdat), tipList[[k]])
				vec[k] <- sum(traitdat[commonsp] != 'Oviparous') / length(commonsp)
				
			}
		}

		if (unlist(selectedTraits2)[i] == 'foragingMode') {
			
			for (k in 1:length(tipList)) {
				
				commonsp <- intersect(names(traitdat), tipList[[k]])
				vec[k] <- sum(traitdat[commonsp] != 'ambush') / length(commonsp)
				
			}
		}

		# if (unlist(selectedTraits2)[i] == 'kinesis') {
			
			# for (k in 1:length(tipList)) {
				
				# commonsp <- intersect(names(traitdat), tipList[[k]])
				# vec[k] <- sum(traitdat[commonsp] != 'akinetic') / length(commonsp)
				
			# }
		# }

		if (unlist(selectedTraits2)[i] == 'prehensionMechanism') {
			
			for (k in 1:length(tipList)) {
				
				commonsp <- intersect(names(traitdat), tipList[[k]])
				vec[k] <- sum(traitdat[commonsp] %in% c('lingual', 'both')) / length(commonsp)
				
			}
		}
				
	} else {

		if (grepl('rate|clads|svl|mass', unlist(selectedTraits2)[i], ignore.case = TRUE)) {
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
ptSizeListRescaled <- lapply(ptSizeList, function(x) rescale(x, to = c(0.5, 2), from = range(unlist(ptSizeList), na.rm = TRUE)))

style <- 'points'

ResizePoints <- FALSE




pdf(paste0('~/Downloads/figS9.pdf'), width = 11, height = 10)

mat <- matrix(c(2,1), nrow = 2, ncol = 1, byrow = TRUE)
layout(mat, heights = c(0.85, 0.15))

par(xpd = NA, oma = c(2,8,3,1))

subtree$root.edge <- 5

qq <- plot.phylo(subtree, show.tip.label = FALSE, edge.width = 0.5, edge.color = 'black', direction = 'upwards', no.margin = TRUE, root.edge = TRUE)

# axisPhylo(side = 2, line = -1, cex.axis = 0.75)
oldest <- max(branching.times(subtree)) + 5
segments(x0 = -2, x1 = -2, y0 = 0, y1 = oldest, xpd = NA, lend = 1, lwd = 0.5)
tx <- seq(from = 0, to = 200, by = 50)
segments(x0 = -3, x1 = -2, y0 = oldest - tx, y1 = oldest - tx, xpd = NA, lend = 1, lwd = 0.5)
text(x = -2.5, y = oldest - tx, labels = tx, pos = 2, cex = 0.5)



# # label major nodes
# toLabel <- c('Serpentes', 'Alethinophidia', 'Colubriformes', 'Gekkota', 'Scincomorpha', 'Anguiformes', 'Iguania', 'Amphisbaenia', 'Lacertoidea')
# adjList <- list(
	# c(-0.1, 1.1),
	# c(-0.1, 1.1),
	# c(-0.1, 1.1),
	# c(-0.1, 1.1),
	# c(1.1, 1.1),
	# c(-0.05, 0.8),
	# c(0, -0.9), # iguania
	# c(-0.05, 1.1), # amphisbaenia
	# c(1.1, 0.9)  # Lacertoidea
# )
# cl <- grep('clade_', colnames(alldat), value = TRUE)
# for (i in 1:length(toLabel)) {
	# nn <- getMRCA(subtree, intersect(subtree$tip.label, alldat[alldat[, paste0('clade_', toLabel[i])] == 1, 'treename']))
	# if (length(nn) > 0) {
		# nodelabels(toLabel[i], node = nn, frame = 'none', cex = 0.9, adj = adjList[[i]])
		# nodelabels(node = nn, pch = 21, cex = 0.5, bg = 'black')
	# }
# }

nodeLabels <- c(A = 'Dibamidae', B = 'Gekkota', C = 'Scincoidea', D = 'Teioidea', E = 'Amphisbaenia', F = 'Lacertoidea', G = 'Anguiformes', H = 'Iguania', I = 'Snakes', J = 'Scolecophidia', K = 'Alethinophidia', L = 'Colubriformes')
points(x = 1, y = 50, pch = 21, bg = 'white', cex = 1, lwd = 0.5)
text(x = 1, y = 50, 'A', cex = 0.3)
for (i in 2:length(nodeLabels)) {
	cc <- paste0('clade_', nodeLabels[i])
	if (nodeLabels[i] == 'Snakes') cc <- 'clade_Serpentes'
	nn <- getMRCA(subtree, intersect(subtree$tip.label, alldat[alldat[, cc] == 1, 'treename']))
	if (length(nn) > 0) {
		nodelabels(node = nn, pch = 21, cex = 1, bg = 'white', lwd = 0.5)
		nodelabels(names(nodeLabels)[i], node = nn, frame = 'none', cex = 0.3)
		
	}
}



# ----------------------
maxPlot <- 1500
int <- seq(300, 1400, length.out = length(unlist(selectedTraits2)))


refW <- diff(int)[1]/2 * 0.85
if ((int[length(int)] + refW) > maxPlot) {
	refW <- maxPlot - int[length(int)]
}

for (i in 2:length(selectedTraits2)) {
	ii <- (which(unlist(selectedTraits2) == selectedTraits2[[i]][1])) : length(int)
	int[ii] <- int[ii] + 30
	segments(x0 = grconvertX(0, from = 'nfc', to ='user'), x1 = Ntip(subtree), y0 = mean(c(int[ii[1]], int[ii[1] - 1])), y1 = mean(c(int[ii[1]], int[ii[1] - 1])), lwd = 0.4)
}

width2 <- 0.4

# axis(2, at = seq(0, 2000, by = 50), xpd = NA, line = -1)
# rect(xleft = 1, ybottom = int - refW, xright = Ntip(subtree), ytop = int + refW, xpd = NA)

rect(xleft = 1, ybottom = int - refW, xright = Ntip(subtree), ytop = int + refW, xpd = NA, col = gray(0.95), border = NA)


colpal <- viridis::plasma
#colpal <- viridis::inferno
# colpal <- colorRampPalette(BAMMtools:::palettes$RdYlBu)
#colpal <- viridis::plasma(5)
#colpal[1] <- 'darkblue'
#colpal <- colorRampPalette(colpal)

colpalDiv <- colorRampPalette(rev(c('#d7191c', '#fdae61','#ffffbf', '#abd9e9', 'darkblue')), alpha = TRUE)

# plot(1:20, rep(1, 20), bg = colpal(20), pch = 21, cex = 5)
# points(1:20, rep(0.8, 20), bg = colpalDiv(20), pch = 21, cex = 5)

segments(
	x0 = rescale(1:100, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	x1 = rescale(1:100, , from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y0 = rescale(rep(0, 99), from = c(0,1), to = c(50, 65)), 
	y1 = rescale(rep(1, 99), from = c(0,1), to = c(50, 65)), 
	col = colpalDiv(100), lend = 1, lwd = 2)
segments(	
	x0 = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	x1 = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y0 = rescale(0, from = c(0,1), to = c(50, 65)), 
	y1 = rescale(1.4, from = c(0,1), to = c(50, 65)), 
	col = 'black', lend = 1, lwd = 0.5)
text(
	x = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(1, from = c(0,1), to = c(50, 65)),
	labels = 'median', pos = 3, cex = 0.75)
text(
	x = rescale(1, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(0.5, from = c(0,1), to = c(50, 65)),
	labels = 'low', pos = 2, cex = 1)
text(
	x = rescale(100, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(0.5, from = c(0,1), to = c(50, 65)),
	labels = 'high', pos = 4, cex = 1)
	
# # scaling legend
# points(x = grconvertX(seq(0.85, 0.92, length.out = 4), from='npc', to = 'user'), y = rep(30, 4), cex = seq(min(unlist(ptSizeListRescaled), na.rm = TRUE), max(unlist(ptSizeListRescaled), na.rm = TRUE), length.out = 4))
# text(x=grconvertX(seq(0.85, 0.92, length.out = 4), from='npc', to = 'user'), y = rep(25, 4), labels = seq(min(unlist(ptSizeList), na.rm = TRUE), max(unlist(ptSizeList), na.rm = TRUE), length.out = 4), pos = 1, cex = 0.8)
# text(x=grconvertX(0.945, from='npc', to = 'user'), y = 25, labels = 'sp.', pos = 1, cex = 0.8)


snakeBounds <- range(which(subtree$tip.label %in% intersect(alldat[alldat$clade_Serpentes == 1, 'treename'], subtree$tip.label)))

colubriformBounds <- range(which(subtree$tip.label %in% intersect(alldat[alldat$clade_Colubriformes == 1, 'treename'], subtree$tip.label)))

sep <- grconvertX(1, from = 'user')

# rect(xleft = snakeBounds[1] - 15, xright = snakeBounds[2] + 0.5, ybottom = 240 - refW, ytop = int[length(int)] + refW+25, xpd = NA, lwd = 0.5, col = adjustcolor(gray(0.9), alpha.f = 0.5), border = NA)

# rect(xleft = colubriformBounds[1] - 10, xright = colubriformBounds[2] + 0.5, ybottom = 240 - refW, ytop = int[length(int)] + refW + 12, xpd = NA, lwd = 0.5, col = adjustcolor(gray(0.9), alpha.f = 0.5), border = NA)

rect(xleft = snakeBounds[1] - 0.5, xright = snakeBounds[2] + 1, ybottom = 240 - 15, ytop = int[length(int)] + refW+35, xpd = NA, lwd = 0.5, lty = 3)

rect(xleft = colubriformBounds[1] - 0.5, xright = colubriformBounds[2] + 1, ybottom = 240 - 15, ytop = int[length(int)] + refW + 12, xpd = NA, lwd = 0.5, lty = 3)


text(x = grconvertX(0.75, from = 'nfc', to ='user'), y = int[length(int)] + refW + 19, 'Colubriformes', pos = 4, xpd = NA, cex = 0.8)
text(x = grconvertX(0.85, from = 'nfc', to ='user'), y = int[length(int)] + refW + 42, 'Serpentes', pos = 4, xpd = NA, cex = 0.8)


refW <- refW * 0.75

#counter <- 1
for (i in 1:length(selectedTraits2)) {

	for (j in 1:length(selectedTraits2[[i]])) {
			
		ij <- which(names(wMeanList) == selectedTraits2[[i]][[j]])
		ref <- int[ij]


		if (grepl('rate|clads', selectedTraits2[[i]][[j]], ignore.case = TRUE)) {

			traitdat <- wMeanDiffList[[ij]]

			refVal <- median(traitdat, na.rm = TRUE)
			refVal <- rescale(refVal, from = range(traitdat, na.rm = TRUE), to = c(ref-refW, ref+refW))
			datRange <- range(traitdat, na.rm = TRUE)

			segments(x0 = 1, x1 = Ntip(subtree), y0 = refVal, y1 = refVal, xpd = NA, lend = 1, lwd = width, lty = 1, col = gray(0.5))
			
			brks <- classIntervals(traitdat, n = 65, style = 'equal')
			cols <- findColours(brks, pal = colpalDiv(length(brks$brks) - 1))

			# lollipops
#			segments(x0 = 1:Ntip(subtree), x1 = 1:Ntip(subtree), y0 = refVal, y1 = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), lwd = 0.15, xpd = NA, lend = 1)
			
		} else {
			
			if (selectedTraits2[[i]][[j]] %in% discrete) {

				traitdat <- wMeanList[[ij]]
				medianDiff <- traitdat
			
			} else {
				
				traitdat <- wMeanList[[ij]]
				medianDiff <- wMeanDiffList[[ij]]

			}
			
#			refVal <- ref-refW
#			datRange <- c(0, max(traitdat, na.rm = TRUE))

			refVal <- median(traitdat, na.rm = TRUE)
			refVal <- rescale(refVal, from = range(traitdat, na.rm = TRUE), to = c(ref-refW, ref+refW))
			datRange <- range(traitdat, na.rm = TRUE)

			segments(x0 = 1, x1 = Ntip(subtree), y0 = refVal, y1 = refVal, xpd = NA, lend = 1, lwd = width, lty = 1, col = gray(0.5))

			brks <- classIntervals(round(medianDiff, 10), n = 65, style = 'equal')
			cols <- findColours(brks, pal = colpalDiv(length(brks$brks) - 1))
			
		}
		
		if (ResizePoints) {
			ptSizes <- ptSizeListRescaled[[ij]]
		} else {
			ptSizes <- 0.8
		}
		
#		# lollipops
#		segments(x0 = 1:Ntip(subtree), x1 = 1:Ntip(subtree), y0 = refVal, y1 = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), lwd = 0.15, xpd = NA, lend = 1)
	
		points(x = 1:Ntip(subtree), y = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), xpd = NA, pch = 21, bg = cols, col = 'black', lwd = 0.15, cex = ptSizes)

		lab <- names(selectedTraits2[[i]])[j]
		if (grepl('chemosensory', lab)) {
			lab <- bquote(psi[chem])
		
		} else if (grepl('vertInnovation', lab)) {
			lab <- 	bquote(psi[vert])

		} else if (grepl('climateRate', lab)) {
			lab <- 	bquote(TR[climate])

		} else if (grepl('vertRate', lab)) {
			lab <- 	bquote(TR[vert])
		}
		
		# if (grepl('rate|clads', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			# lab <- bquote(beta[.(names(selectedTraits[[i]])[j])])
		# } else if (grepl('ancDist|chemo', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			# lab <- bquote(Phi[.(names(selectedTraits[[i]])[j])])
		# } else {
			# lab <- names(selectedTraits[[i]])[j]
		# }
		text(x = grconvertX(0.03, from = 'nfc', to ='user'), y = ref, labels = lab, xpd = NA, pos = 2, srt = 0, cex = 0.75)
		
		if (j == 1 & i != 1) {
			text(x = grconvertX(-0.05, from = 'nfc', to ='user'), y = mean(int[unlist(selectedTraits2) %in% selectedTraits2[[i]]]), labels = names(selectedTraits2)[i], xpd = NA, pos = 3, srt = 90)
		}
	}	
}

# secondary, type 1
segments(x0 = grconvertX(-0.1, from = 'nfc', to ='user'), x1 = grconvertX(-0.1, from = 'nfc', to ='user'), y0 = grconvertY(0.94, from = 'ndc'), y1 = grconvertY(0.78, from = 'ndc'))
text(x = grconvertX(-0.12, from = 'nfc', to ='user'), y = grconvertY(0.86, from = 'ndc'), labels = 'Type 1 traits', xpd = NA, pos = 3, srt = 90, font = 2)

# secondary type 2
segments(x0 = grconvertX(-0.1, from = 'nfc', to ='user'), x1 = grconvertX(-0.1, from = 'nfc', to ='user'), y0 = grconvertY(0.75, from = 'ndc'), y1 = grconvertY(0.39, from = 'ndc'))
text(x = grconvertX(-0.12, from = 'nfc', to ='user'), y = grconvertY(0.57, from = 'ndc'), labels = 'Type 2 traits', xpd = NA, pos = 3, srt = 90, font = 2)

# biogeographic
segments(x0 = grconvertX(-0.1, from = 'nfc', to ='user'), x1 = grconvertX(-0.1, from = 'nfc', to ='user'), y0 = grconvertY(0.36, from = 'ndc'), y1 = grconvertY(0.2, from = 'ndc'))
text(x = grconvertX(-0.12, from = 'nfc', to ='user'), y = grconvertY(0.28, from = 'ndc'), labels = 'biogeographic', xpd = NA, pos = 3, srt = 90, font = 2)


#text(x = grconvertX(-0.12, from = 'nfc', to ='user'), y = grconvertY(0.5, from = 'ndc'), labels = 'primary traits', xpd = NA, pos = 3, srt = 90, font = 2)

dev.off()




