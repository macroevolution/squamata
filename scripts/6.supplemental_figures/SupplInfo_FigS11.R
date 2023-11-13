
# histograms of rate differences

# elongation index and rate
# vertebral counts and rate
# skull shape and rate

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

require(ape)

source('./scripts/6.supplemental_figures/violin.R')

alldat <- read.csv('alldat.csv')

tr <- read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")

snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']

f3 <- './data/3.trait_preparation_diversification/2D_Adult_Skull_Coord_treenames.csv'
skullCoords 		<- read.csv(f3)
rownames(skullCoords) <- skullCoords$treename
skullCoords <- skullCoords[, grep('ProcCoord', colnames(skullCoords))]
# skullCoords <- skullCoords[intersect(rownames(skullCoords), tr$tip.label),]
head(skullCoords)
skullpca <- prcomp(skullCoords)
summary(skullpca)


snakeInd <- which((alldat$treename %in% snakes) == TRUE)
lizardInd <- which((alldat$treename %in% lizards) == TRUE)

pdf('~/Downloads/figS11.pdf', width = 10, height = 5)

par(mfcol = c(2, 3))

par(mar = c(5,5,1,1), oma = c(0.5, 0.5, 1, 0.5), lwd = 0.5)

labelSize <- 0.7
letterHeight <- -0.5

# elongation
hist(log(alldat[lizardInd, 'elongationIndex']), col = adjustcolor('slateblue1', alpha.f = 0.75), axes = FALSE, xlim = c(4.7, 8), xlab = '', ylab = '', breaks = 'fd', main = NA)
hist(log(alldat[snakeInd, 'elongationIndex']), col = adjustcolor('coral2', alpha.f = 0.75), breaks = 'fd', add = TRUE)
axis(1, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
axis(2, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
mtext('log elongation index', side = 1, line = 2, cex = labelSize)
mtext('species count', side = 2, line = 2, cex = labelSize)
box(which = "plot", bty = "l")
mtext('A', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

mtext('elongation index', side = 3)

legend('right', legend = c('snakes', 'lizards'), fill = adjustcolor(c('coral2', 'slateblue1'), alpha.f = 0.75), bty = 'n')

hist(log(alldat[lizardInd, 'elongationRate']), col = adjustcolor('slateblue1', alpha.f = 0.75), axes = FALSE, xlim = c(-16, 2), breaks = 'fd', xlab = '', ylab = '', main = NA)
hist(log(alldat[snakeInd, 'elongationRate']), col = adjustcolor('coral2', alpha.f = 0.75), breaks = 'fd', add = TRUE)
axis(1, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
axis(2, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
mtext('log elongation rate', side = 1, line = 2, cex = labelSize)
mtext('species count', side = 2, line = 2, cex = labelSize)
box(which = "plot", bty = "l")
mtext('D', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)


# vertebral counts
hist(log(alldat[lizardInd, 'numberPresacralVert']), col = adjustcolor('slateblue1', alpha.f = 0.75), axes = FALSE, xlim = c(2.5, 6), breaks = 'fd', xlab = '', ylab = '', main = NA)
hist(log(alldat[snakeInd, 'numberPresacralVert']), col = adjustcolor('coral2', alpha.f = 0.75), breaks = 'fd', add = TRUE)
axis(1, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
axis(2, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
mtext('log presacral vertebral count', side = 1, line = 2, cex = labelSize)
mtext('species count', side = 2, line = 2, cex = labelSize)
box(which = "plot", bty = "l")
mtext('B', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

mtext('vertebral counts', side = 3)

hist(log(alldat[lizardInd, 'vertebralRate']), col = adjustcolor('slateblue1', alpha.f = 0.75), axes = FALSE, xlim = c(-7.5, 9.5), ylim = c(0, 150), breaks = 'fd', xlab = '', ylab = '', main = NA)
hist(log(alldat[snakeInd, 'vertebralRate']), col = adjustcolor('coral2', alpha.f = 0.75), breaks = 'fd', add = TRUE)
axis(1, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
axis(2, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
mtext('log presacral vertebral count rate', side = 1, line = 2, cex = labelSize)
mtext('species count', side = 2, line = 2, cex = labelSize)
box(which = "plot", bty = "l")
mtext('E', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)


plot.new()
plot.window(xlim = range(skullpca$x[,1]), ylim = range(skullpca$x[,2]))
axis(1, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
axis(2, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
box(which = "plot", bty = "l")
mtext(paste0('skull shape PC 1 (', round(summary(skullpca)$importance[2,1]*100,1), '%)'), side = 1, line = 2, cex = labelSize)
mtext(paste0('skull shape PC 2 (', round(summary(skullpca)$importance[2,2]*100,1), '%)'), side = 2, line = 2, cex = labelSize)
cols <- rep('slateblue1', nrow(skullpca$x))
cols[rownames(skullpca$x) %in% snakes] <- 'coral2'
points(skullpca$x[,1], skullpca$x[,2], bg = cols, lwd = 0.75, pch = 21)
mtext('C', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)


mtext('skull shape', side = 3)

# breaks 20, 15
hist(log(alldat[lizardInd, 'skullRate']), col = adjustcolor('slateblue1', alpha.f = 0.75), axes = FALSE, xlim = c(-10, -5), breaks = 'fd', xlab = '', ylab = '', main = NA)
hist(log(alldat[snakeInd, 'skullRate']), col = adjustcolor('coral2', alpha.f = 0.75), breaks = 'fd', add = TRUE)
axis(1, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
axis(2, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
mtext('log skull shape rate', side = 1, line = 2, cex = labelSize)
mtext('species count', side = 2, line = 2, cex = labelSize)
box(which = "plot", bty = "l")
mtext('F', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

dev.off()

