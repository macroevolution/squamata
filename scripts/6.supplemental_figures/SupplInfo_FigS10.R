## Figure S9: traits, rates, net innovation but for all species

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


## Supplemental figure showing trait data but for all species
selectedTraits <- 	list(
						"Rates" = list(
							'CLaDS' = 'meanImputedCLADS',
							'climate' = 'climateNicheRate',
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
						),
						"Raw traits" = list(
							'climate\nPC1' = 'climPC1',
							'diet\nPC1' = 'dietPCA',
							'skull shape\nPC1' = 'skullPCA',
							'log elongation\nindex' = 'elongationIndex',
							'log mass' = 'mass',
							'log # presacral\nvertebrae' = 'numberPresacralVert'
						)
					)

discrete <- c('parity', 'numberDigits', 'numberLimbs', 'kinesis', 'prehensionMechanism', 'foragingMode')


subtree <- tr
	
f3 <- './data/3.trait_preparation_diversification/2D_Adult_Skull_Coord_treenames.csv'
skullCoords 		<- read.csv(f3)
rownames(skullCoords) <- skullCoords$treename
skullCoords <- skullCoords[, grep('ProcCoord', colnames(skullCoords))]
skullCoords <- skullCoords[intersect(rownames(skullCoords), subtree$tip.label),]
head(skullCoords)
skullpca <- prcomp(skullCoords)
summary(skullpca)

# diet proportions
dietProportions <- read.csv('./data/3.trait_preparation_diversification/dietModeling/diet-proportions-phylo.csv', row.names = 1)
dietProportions <- data.matrix(dietProportions)
# dietProportions <- log(sweep(dietProportions[,-31], 1, dietProportions[,31], "/")) # alr
dietProportions <- sweep(log(dietProportions), 1, rowMeans(log(dietProportions)))
dietpca = prcomp(dietProportions)
summary(dietpca)



library(RColorBrewer)


png('~/Downloads/figS10.png', width = 10, height = 8, units = 'in', res = 600)

mat <- matrix(c(2,1), nrow = 2, ncol = 1, byrow = TRUE)
layout(mat, heights = c(0.90, 0.10))

par(oma = c(1,5,1,1), xpd = NA)


qq <- plot.phylo(subtree, show.tip.label = FALSE, edge.width = width, cex = 0.12, edge.color = 'black', direction = 'upwards', no.margin = TRUE)

# axisPhylo(side = 2, line = -1, cex.axis = 0.75)
oldest <- max(branching.times(subtree))
segments(x0 = -20, x1 = -20, y0 = 0, y1 = oldest, xpd = NA, lend = 1, lwd = 0.5)
tx <- seq(from = 0, to = 200, by = 50)
segments(x0 = -40, x1 = -20, y0 = oldest - tx, y1 = oldest - tx, xpd = NA, lend = 1, lwd = 0.5)
text(x = -2.5, y = oldest - tx, labels = tx, pos = 2, cex = 0.4)



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
maxPlot <- 2200
int <- seq(280, 2000, length.out = length(unlist(selectedTraits)))


refW <- diff(int)[1]/2 * 0.85
if ((int[length(int)] + refW) > maxPlot) {
	refW <- maxPlot - int[length(int)]
}

for (i in 2:length(selectedTraits)) {
	ii <- (which(unlist(selectedTraits) == selectedTraits[[i]][1])) : length(int)
	int[ii] <- int[ii] + 30
	segments(x0 = grconvertX(-0.02, from = 'nfc', to ='user'), x1 = Ntip(subtree), y0 = mean(c(int[ii[1]], int[ii[1] - 1])), y1 = mean(c(int[ii[1]], int[ii[1] - 1])), lwd = 0.4, xpd = NA)
}

width2 <- 0.4

# axis(2, at = seq(0, 2200, by = 50), xpd = NA, line = -1)
# rect(xleft = 1, ybottom = int - refW, xright = Ntip(subtree), ytop = int + refW, xpd = NA)

rect(xleft = 1, ybottom = int - refW, xright = Ntip(subtree), ytop = int + refW, xpd = NA, col = gray(0.95), border = NA)

snakeBounds <- range(which(subtree$tip.label %in% intersect(alldat[alldat$clade_Serpentes == 1, 'treename'], subtree$tip.label)))

sep <- grconvertX(1, from = 'user')

rect(xleft = snakeBounds[1] - 0.5, xright = snakeBounds[2] + 1, ybottom = 240 - 25, ytop = int[length(int)] + refW+25, xpd = NA, lwd = 0.5, lty = 3)

text(x = grconvertX(0.85, from = 'nfc', to ='user'), y = int[length(int)] + refW + 38, 'Serpentes', pos = 4, xpd = NA, cex = 0.7)


colpal <- viridis::plasma

colpalDiv <- colorRampPalette(rev(c('#d7191c', '#fdae61','#ffffbf', '#abd9e9', 'darkblue')), alpha = TRUE)

# plot(1:20, rep(1, 20), bg = colpal(20), pch = 21, cex = 5)
# points(1:20, rep(0.8, 20), bg = colpalDiv(20), pch = 21, cex = 5)

segments(
	x0 = rescale(1:100, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	x1 = rescale(1:100, , from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y0 = rescale(rep(0, 99), from = c(0,1), to = c(70, 85)), 
	y1 = rescale(rep(1, 99), from = c(0,1), to = c(70, 85)), 
	col = colpalDiv(100), lend = 1, lwd = 2)
segments(	
	x0 = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	x1 = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y0 = rescale(0, from = c(0,1), to = c(70, 85)), 
	y1 = rescale(1.4, from = c(0,1), to = c(70, 85)), 
	col = 'black', lend = 1, lwd = 0.5)
text(
	x = rescale(50, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(1, from = c(0,1), to = c(70, 85)),
	labels = 'median', pos = 3, cex = 0.65)
text(
	x = rescale(1, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(0.5, from = c(0,1), to = c(70, 85)),
	labels = 'low', pos = 2, cex = 0.75)
text(
	x = rescale(100, from = c(1,100), to = grconvertX(c(0.85, 0.92), from = 'nfc', to ='user')), 
	y = rescale(0.5, from = c(0,1), to = c(70, 85)),
	labels = 'high', pos = 4, cex = 0.75)


refW <- refW * 0.75

for (i in 1:length(selectedTraits)) {
	
	for (j in 1:length(selectedTraits[[i]])) {
			
		ij <- which(unlist(selectedTraits) == selectedTraits[[i]][[j]])
		ref <- int[ij]
		
		if (selectedTraits[[i]][[j]] == 'skullPCA') {
			
			traitdat <- rep(NA, Ntip(subtree))
			names(traitdat) <- subtree$tip.label
			traitdat[match(rownames(skullpca$x), names(traitdat))] <- skullpca$x[,1]
	
		} else if (selectedTraits[[i]][[j]] == 'dietPCA') {
			
			traitdat <- rep(NA, Ntip(subtree))
			names(traitdat) <- subtree$tip.label
			traitdat[match(rownames(dietpca$x), names(traitdat))] <- dietpca$x[,1]
	
	
		} else {
	
			traitdat <- setNames(alldat[, selectedTraits[[i]][[j]]], alldat$treename)[subtree$tip.label]
			# range(traitdat, na.rm = TRUE)
	
		}
		
		if ((selectedTraits[[i]][[j]] %in% c('meanImputedCLADS', 'completeSVL', 'mass', 'elongationIndex', 'dietBreadth', 'elongationIndex', 'rangeSize', 'numberPresacralVert')) | grepl('rate', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			
			if (min(traitdat, na.rm = TRUE) == 0) {
				traitdat <- traitdat + 0.01
			}
			traitdat <- log(traitdat)
		}
		
		if (selectedTraits[[i]][[j]] %in% discrete) {

			traitdat <- as.numeric(traitdat)
			brks <- classIntervals(traitdat, n = 65, style = 'equal')
			cols <- findColours(brks, pal = colorRampPalette(c('#1b9e77', '#d95f02', '#7570b3'))(length(brks$brks)-1))

			points(x = 1:Ntip(subtree), y = rescale(traitdat, to = c(ref-refW, ref+refW)), xpd = NA, pch = 21, bg = cols, col = 'black', lwd = 0.05, cex = 0.1)
				
			
		} else {
		
			medianTraitDat <- median(traitdat, na.rm = TRUE)
			traitDiffs <- traitdat - medianTraitDat
				
			if (grepl('rate|clads', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
		
				refVal <- 0
				refVal <- rescale(refVal, from = range(traitDiffs, na.rm = TRUE), to = c(ref-refW, ref+refW))
				datRange <- range(traitDiffs, na.rm = TRUE)
	
				segments(x0 = 1, x1 = Ntip(subtree), y0 = refVal, y1 = refVal, xpd = NA, lend = 1, lwd = width, lty = 1, col = gray(0.5))
				
				brks <- classIntervals(traitDiffs, n = 65, style = 'equal')
				cols <- findColours(brks, pal = colpalDiv(length(brks$brks) - 1))
				
				points(x = 1:Ntip(subtree), y = rescale(traitDiffs, from = datRange, to = c(ref-refW, ref+refW)), xpd = NA, pch = 21, bg = cols, col = 'black', lwd = 0.05, cex = 0.1)
			
			} else {
				
				refVal <- medianTraitDat
				if (selectedTraits[[i]][[j]] == 'centroidLat') refVal <- 0
				refVal <- rescale(refVal, from = range(traitdat, na.rm = TRUE), to = c(ref-refW, ref+refW))
				datRange <- range(traitdat, na.rm = TRUE)
	
				segments(x0 = 1, x1 = Ntip(subtree), y0 = refVal, y1 = refVal, xpd = NA, lend = 1, lwd = width, lty = 1, col = gray(0.5))
	
				brks <- classIntervals(traitDiffs, n = 65, style = 'equal')
				cols <- findColours(brks, pal = colpalDiv(length(brks$brks) - 1))
					
				points(x = 1:Ntip(subtree), y = rescale(traitdat, from = datRange, to = c(ref-refW, ref+refW)), xpd = NA, pch = 21, bg = cols, col = 'black', lwd = 0.05, cex = 0.1)
				
				
			}
		}
		
		if (grepl('ancDist|chemo', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			lab <- bquote(psi[.(names(selectedTraits[[i]])[j])])
		} else if (grepl('rate', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			lab <- bquote(TR[.(names(selectedTraits[[i]])[j])])
		} else if (grepl('clads|bamm', selectedTraits[[i]][[j]], ignore.case = TRUE)) {
			lab <- bquote(lambda[.(names(selectedTraits[[i]])[j])])
		} else {
			lab <- names(selectedTraits[[i]])[j]
		}
	
			
		text(x = grconvertX(0.04, from = 'nfc', to ='user'), y = ref, labels = lab, xpd = NA, pos = 2, srt = 0, cex = 0.5)
		
		if (j == 1) {
			text(x = grconvertX(-0.04, from = 'nfc', to ='user'), y = mean(int[unlist(selectedTraits) %in% selectedTraits[[i]]]), labels = names(selectedTraits)[i], xpd = NA, pos = 3, srt = 90, cex = 0.7)
		}
	}	
}

dev.off()

