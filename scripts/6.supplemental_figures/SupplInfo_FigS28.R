# latitude, richness and speciation rates


setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

require(ape)
require(scales)
require(classInt)


alldat <- read.csv('./data/alldat.csv')

treefile <- './data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre'

tr <- read.tree(treefile)
tr <- ladderize(tr, right = FALSE)

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
tr <- rotate(tr, nn)

tr <- read.tree(text = write.tree(tr)) # otherwise, node rotation doesn't change structure



snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']

# explore: remove Liolaemidae
# lizards <- setdiff(lizards, alldat[alldat$family == 'Liolaemidae', 'treename'])

table(is.na(alldat$centroidLat))

alldat <- alldat[!is.na(alldat$centroidLat), ]

latBins <- seq(-90, 90, by = 1)

#latList <- vector('list', length(latBins))
#names(latList) <- as.character(latBins)
resmat <- matrix(nrow = length(latBins), ncol = 6)
colnames(resmat) <- c('allRichness', 'lizardRichness', 'snakeRichness', 'allRates', 'lizardRates', 'snakeRates')
resmat <- as.data.frame(resmat)
for (i in 1:length(latBins)) {
	
	# centroid
	#latList[[i]] <- alldat[alldat$centroidLat > latBins[i] & alldat$centroidLat <= latBins[i+1], 'treename']
	
	# using latitudinal range
	allsp <- alldat[alldat$maxLat >= latBins[i] & alldat$minLat <= latBins[i+1], 'treename']
	resmat[i, 'allRichness'] <- length(allsp)	
	resmat[i, 'lizardRichness'] <- length(intersect(allsp, lizards))
	resmat[i, 'snakeRichness'] <- length(intersect(allsp, snakes))

	# mean imputed CLaDS
	resmat[i, 'allRates'] <- mean(alldat[alldat$treename %in% allsp, 'meanImputedCLADS'], na.rm = TRUE)
	resmat[i, 'lizardRates'] <- mean(alldat[alldat$treename %in% intersect(allsp, lizards), 'meanImputedCLADS'], na.rm = TRUE)
	resmat[i, 'snakeRates'] <- mean(alldat[alldat$treename %in% intersect(allsp, snakes), 'meanImputedCLADS'], na.rm = TRUE)
	
}

# binRate <- sapply(latList, function(x) mean(alldat[alldat$treename %in% x, 'meanImputedCLADS'], na.rm = TRUE))

snakeCol <- adjustcolor('coral2', alpha.f = 0.5)
lizardCol <- adjustcolor('slateblue1', alpha.f = 0.5)

axisLabSize <- 0.75
ptSize <- 0.65

plotRanges <- TRUE
plotPoints <- FALSE



#####################################

# VERTICAL

tr <- read.tree(treefile)
tr <- ladderize(tr, right = TRUE)

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
tr <- rotate(tr, nn)

tr <- read.tree(text = write.tree(tr)) # otherwise, node rotation doesn't change structure

axisLabSize <- 0.75
ptSize <- 0.75


pdf('~/Downloads/figS28.pdf', width = 8, height = 9)

mat <- matrix(c(4,5,1,4,5,2,4,5,3), nrow = 3, ncol = 3, byrow = TRUE)
layout(mat, widths = c(0.2, 0.45, 0.45))

par(xpd = NA, oma = c(1,1,1,1), mar = c(4, 4, 2, 2))


# sp richness ~ latitude
plot.new()
plot.window(xlim = range(latBins), ylim = c(0, 800))
axis(1, at = seq(-80, 80, by = 40), lwd = 0, lwd.ticks = 1, cex.axis = 0.75)
axis(2, lwd = 0, lwd.ticks = 1, las = 2, cex.axis = 0.75)
box(which = "plot", bty = "l")
#points(latBins, resmat$allRichness)
points(latBins, resmat$lizardRichness, bg = lizardCol, lwd = 0.5, pch = 21, cex = ptSize)
points(latBins, resmat$snakeRichness, bg = snakeCol, lwd = 0.5, pch = 21, cex = ptSize)
mtext('latitude', side = 1, line = 2.5, cex = axisLabSize)
mtext('sp richness', side = 2, line = 2.5, cex = axisLabSize)
text(x = grconvertX(-0.1, from = 'npc', to = 'user'), y = grconvertY(1.1, from = 'npc', to = 'user'), 'B', xpd = NA, cex = 1.2, font = 2)

legend('topright', pt.bg = c(lizardCol, snakeCol), legend = c('lizards', 'snakes'), pt.lwd = 0.5, pch = 21, bty = 'n')

# speciation rate ~ latitude
ind <- which(resmat$allRichness >= 5)
plot.new()
plot.window(xlim = range(latBins), ylim = range(log(resmat[ind, c('lizardRates', 'snakeRates')]), na.rm = TRUE))
axis(1, at = seq(-80, 80, by = 40), lwd = 0, lwd.ticks = 1, cex.axis = 0.75)
axis(2, lwd = 0, lwd.ticks = 1, cex.axis = 0.75, las = 2)
box(which = "plot", bty = "l")
# points(latBins, resmat$allRates)
points(latBins[ind], log(resmat$lizardRates[ind]), bg = lizardCol, lwd = 0.5, pch = 21, cex = ptSize)
points(latBins[ind], log(resmat$snakeRates[ind]), bg = snakeCol, lwd = 0.5, pch = 21, cex = ptSize)
mtext('latitude', side = 1, line = 2.5, cex = axisLabSize)
mtext('log speciation rate', side = 2, line = 2.5, cex = axisLabSize)
text(x = grconvertX(-0.1, from = 'npc', to = 'user'), y = grconvertY(1.1, from = 'npc', to = 'user'), 'C', xpd = NA, cex = 1.2, font = 2)

# richness ~ speciation rate
plot.new()
plot.window(xlim = log(c(1, max(resmat[ind, c('lizardRichness', 'snakeRichness')]))), ylim = range(log(resmat[ind, c('lizardRates', 'snakeRates')]), na.rm = TRUE), )
axis(1, lwd = 0, lwd.ticks = 1, cex.axis = 0.75)
axis(2, lwd = 0, lwd.ticks = 1, cex.axis = 0.75, las = 2)
box(which = "plot", bty = "l")
# points(binRate, resmat$allRichness)
points(log(resmat$lizardRichness)[ind], log(resmat$lizardRates)[ind], bg = lizardCol, lwd = 0.5, pch = 21, cex = ptSize)
points(log(resmat$snakeRichness)[ind], log(resmat$snakeRates)[ind], bg = snakeCol, lwd = 0.5, pch = 21, cex = ptSize)

mtext('log sp richness', side = 1, line = 2.5, cex = axisLabSize)
mtext('log sspeciation rate', side = 2, line = 2.5, cex = axisLabSize)
text(x = grconvertX(-0.1, from = 'npc', to = 'user'), y = grconvertY(1.1, from = 'npc', to = 'user'), 'D', xpd = NA, cex = 1.2, font = 2)


par(mar = c(1, 1, 1, 1))

tr$root.edge <- 5

qq <- plot.phylo(tr, show.tip.label = FALSE, edge.width = 0.1, edge.color = 'black', direction = 'right', no.margin = TRUE, root.edge = TRUE)

# axisPhylo(side = 2, line = -2.5, cex.axis = 0.75)

# # label major nodes
# toLabel <- c('Serpentes', 'Macrostomata', 'Colubriformes', 'Gekkota', 'Scincoidea', 'Anguimorpha', 'Iguania', 'Amphisbaenia', 'Lacertoidea')
# cl <- grep('clade_', colnames(alldat), value = TRUE)
# for (i in 1:length(toLabel)) {
	# nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, paste0('clade_', toLabel[i])] == 1, 'treename']))
	# if (length(nn) > 0) {
		# nodelabels(toLabel[i], node = nn, frame = 'none', cex = 1, adj = c(-0.1, 1.1))
	# }
# }


zeroLat <- 500
width <- 250

colpal <- colorRampPalette(rev(c('#d7191c', '#fdae61','#ffffbf', '#abd9e9', 'darkblue')), alpha = TRUE)
colpal <- viridis::plasma

latRates <- setNames(alldat[, 'meanImputedCLADS'], alldat$treename)
latRates <- log(latRates[tr$tip.label])
brks <- classIntervals(latRates, n = 65, style = 'equal')
cols <- findColours(brks, pal = colpal(length(brks$brks) - 1))

if (!plotRanges) {
	# axis(1, at = seq(0, 2000, by = 50), xpd = NA, line = -2)
	rect(ybottom = 1, xleft = zeroLat - width, ytop = Ntip(tr), xright = zeroLat + width, xpd = NA)
}

# latitude lines
tropLats <- rescale(c(-23.5, 23.5), from = c(-90,90), to = c(zeroLat - width, zeroLat + width))
segments(y0 = 1, y1 = Ntip(tr), x0 = tropLats, x1 = tropLats, lwd = 0.5, lty = 2)
segments(y0 = 1, y1 = Ntip(tr), x0 = zeroLat, x1 = zeroLat, lwd = 0.5)

# axis bar
segments(y0 = -75, y1 = -75, x0 = rescale(-90, from = c(-90,90), to = c(zeroLat - width, zeroLat + width)), x1 = rescale(90, from = c(-90,90), to = c(zeroLat - width, zeroLat + width)), xpd = NA, lend = 1, lwd = 0.5)

segments(y0 = -110, y1 = -75, x0 = rescale(seq(-90, 90, by = 30), from = c(-90,90), to = c(zeroLat - width, zeroLat + width)), x1 = rescale(seq(-90, 90, by = 30), from = c(-90,90), to = c(zeroLat - width, zeroLat + width)), xpd = NA, lend = 1, lwd = 0.5)

text(y = -100, x = rescale(seq(-90, 90, by = 30), from = c(-90,90), to = c(zeroLat - width, zeroLat + width)), as.character(seq(-90, 90, by = 30)), pos = 1, xpd = NA)

text(y = -180, x = zeroLat, 'latitude', pos = 1, xpd = NA, cex = 1)

# snake box
snakeBounds <- range(which(tr$tip.label %in% intersect(snakes, tr$tip.label)))

sep <- grconvertY(1, from = 'user')

rect(ybottom = snakeBounds[1] - 1, ytop = snakeBounds[2] + 1, xleft = rescale(-70, from = c(-90, 90), to = c(zeroLat - width, zeroLat + width)), xright = rescale(70, from = c(-90, 90), to = c(zeroLat - width, zeroLat + width)), xpd = NA, lwd = 0.5, lty = 3)

text(y = grconvertY(0.3, from = 'nfc', to ='user'), x = rescale(70, from = c(-90, 90), to = c(zeroLat - width, zeroLat + width)), 'Serpentes', pos = 4, xpd = NA, cex = 0.8, srt = 270)

text(x = grconvertX(1.25, from = 'npc', to = 'user'), y = grconvertY(0.995, from = 'npc', to = 'user'), 'A', xpd = NA, cex = 1.2, font = 2)

if (plotRanges) {
	for (i in 1:Ntip(tr)) {
		
		if (tr$tip.label[i] %in% alldat$treename) {
			
			spdat <- alldat[alldat$treename == tr$tip.label[i], ]
			segments(y0 = i, y1 = i, x0 = rescale(spdat[1, 'minLat'], from = c(-90,90), to = c(zeroLat - width, zeroLat + width)), x1 = rescale(spdat[1, 'maxLat'], from = c(-90,90), to = c(zeroLat - width, zeroLat + width)), xpd = NA, lend = 1, lwd = 0.25, col = cols[i])
		}
	}
}

segments(
	y0 = rescale(1:100, from = c(1,100), to = grconvertY(c(0.45, 0.55), from = 'nfc', to ='user')), 
	y1 = rescale(1:100, , from = c(1,100), to = grconvertY(c(0.45, 0.55), from = 'nfc', to ='user')), 
	x0 = rescale(-90, from = c(-90, 90), to = c(zeroLat - width, zeroLat + width)), 
	x1 = rescale(-80, from = c(-90, 90), to = c(zeroLat - width, zeroLat + width)), 
	col = colpal(100), lend = 1, lwd = 2)
text(
	y = rescale(50, from = c(1,100), to = grconvertY(c(0.45, 0.55), from = 'nfc', to ='user')), 
	x = rescale(-77, from = c(-90, 90), to = c(zeroLat - width, zeroLat + width)),
	labels = 'log speciation rates', pos = 3, cex = 0.8, srt = 270)
text(
	y = rescale(1, from = c(1,100), to = grconvertY(c(0.45, 0.55), from = 'nfc', to ='user')), 
	x = rescale(-85, from = c(-90, 90), to = c(zeroLat - width, zeroLat + width)),
	labels = 'slow', pos = 1, cex = 0.8)
text(
	y = rescale(100, from = c(1,100), to = grconvertY(c(0.45, 0.55), from = 'nfc', to ='user')), 
	x = rescale(-85, from = c(-90, 90), to = c(zeroLat - width, zeroLat + width)),
	labels = 'fast', pos = 3, cex = 0.8)
	




dev.off()


