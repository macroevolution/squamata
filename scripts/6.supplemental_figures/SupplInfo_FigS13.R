# body mass variation

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

require(ape)

source('./scripts/6.supplemental_figures/violin.R')

alldat <- read.csv('./data/alldat.csv')

tr <- read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")

snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']

snakeInd <- which((alldat$treename %in% snakes) == TRUE)
lizardInd <- which((alldat$treename %in% lizards) == TRUE)



pdf('~/Downloads/figS13.pdf', width = 12, height = 5)

par(mfrow = c(1,2), lwd = 0.5, c(3,3,1,0.5), oma = c(0, 0, 0, 0))

# breaks 30, 30
hist(log(alldat[lizardInd, 'mass'], base = 10), col = adjustcolor('slateblue1', alpha.f = 0.75), xlim = c(-1, 6), axes = FALSE, breaks = 'fd', xlab = '', ylab = '', main = NA)
hist(log(alldat[snakeInd, 'mass'], base = 10), col = adjustcolor('coral2', alpha.f = 0.75), breaks = 'fd', add = TRUE)
axis(1, at = log(c(1, 10, 100, 1000, 10000, 100000), base = 10), labels = c(1, 10, 100, 1000, 10000, 100000), lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
axis(2, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
# mtext('log mass (g)', side = 1, line = 2, cex = 0.75)
mtext(expression(log[10]~mass~(g)), side = 1, line = 2, cex = 0.75)
mtext('species count', side = 2, line = 2, cex = 0.75)
box(which = "plot", bty = "l")
mtext('A', adj = -0.1, padj = -1, font = 2, xpd = NA, cex = 1.2)


mtext('body mass', side = 3)


hist(log(alldat[lizardInd, 'massRate']), col = adjustcolor('slateblue1', alpha.f = 0.75), xlim = c(-7.5, 6.5), axes = FALSE, breaks = 'fd', xlab = '', ylab = '', main = NA)
hist(log(alldat[snakeInd, 'massRate']), col = adjustcolor('coral2', alpha.f = 0.75), breaks = 'fd', add = TRUE)
axis(1, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
axis(2, lwd = 0, lwd.ticks = 0.5, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0), lend = 1)
mtext('log mass tip rate', side = 1, line = 2, cex = 0.75)
mtext('species count', side = 2, line = 2, cex = 0.75)
box(which = "plot", bty = "l")
mtext('B', adj = -0.1, padj = -1, font = 2, xpd = NA, cex = 1.2)

mtext('body mass rate', side = 3)


dev.off()
