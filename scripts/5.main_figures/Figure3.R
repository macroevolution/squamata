setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

library(ape)
library(phytools)
library(TeachingDemos)
library(scales)

source('./scripts/6.supplemental_figures/violin.R')


alldat <- read.csv('./data/alldat.csv')

tree <- read.tree('./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre')
tree <- read.tree(text = write.tree(ladderize(tree)))

snakes <- alldat[which(alldat[, 'clade_Serpentes'] == 1), 'treename']
snakes <- snakes[snakes %in% tree$tip.label]
lizards <- alldat[which(alldat[, 'clade_Serpentes'] == 0), 'treename']
lizards <- lizards[lizards %in% tree$tip.label]


# diet breadth
breadth <- setNames(alldat$dietBreadth, alldat$treename)
breadth <- breadth[!is.na(breadth)]


# diet proportions
dietProportions <- read.csv('./data/3.trait_preparation_diversification/dietModeling/diet-proportions-phylo.csv', row.names = 1)
dietProportions <- data.matrix(dietProportions)
# dietProportions <- log(sweep(dietProportions[,-31], 1, dietProportions[,31], "/")) # alr
dietProportionsStan <- sweep(log(dietProportions), 1, rowMeans(log(dietProportions))) # clr

pca = prcomp(dietProportionsStan)

ptCols <- rep('slateblue1', nrow(pca$x))
ptCols[rownames(pca$x) %in% snakes] <- 'coral2'
table(ptCols)

ptSizes <- scales::rescale(breadth[rownames(pca$x)], to = c(0.5, 2.5))

labLineCol <- gray(0.75)


par(xpd = NA, mar = c(4,4,1,0))

plot.new()
plot.window(xlim = range(pca$x[,1]), ylim = range(pca$x[,2]))

axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")
mtext(paste0('diet PC 1 (', round(summary(pca)$importance[2,1]*100,1), '%)'), side = 1, line = 2.5)
mtext(paste0('diet PC 2 (', round(summary(pca)$importance[2,2]*100,1), '%)'), side = 2, line = 2.5)

blindsnakes <- alldat[alldat$family %in% c('Typhlopidae', 'Leptotyphlopidae'), 'treename']
blindsnakes <- intersect(blindsnakes, rownames(pca$x))
ptLab <- c(-10, -2)
segments(x0 = pca$x[blindsnakes, 1], x1 = rep(ptLab[1], length(blindsnakes)), y0 = pca$x[blindsnakes, 2], y1 = rep(ptLab[2], length(blindsnakes)), lend = 1, lwd = 0.5, col = labLineCol)
text(ptLab[1], ptLab[2], 'blind snakes', pos = 2, font = 2, col = 'coral2')

heloderma <- alldat[alldat$family %in% c('Helodermatidae'), 'treename']
heloderma <- intersect(heloderma, rownames(pca$x))
ptLab <- c(-2, -1)
segments(x0 = pca$x[heloderma, 1], x1 = rep(ptLab[1], length(heloderma)), y0 = pca$x[heloderma, 2], y1 = rep(ptLab[2], length(heloderma)), lend = 1, lwd = 0.5, col = labLineCol)
text(ptLab[1], ptLab[2], 'Heloderma', pos = 1, font = 4, col = 'slateblue1')

varanus <- alldat[alldat$family %in% c('Varanidae'), 'treename']
varanus <- intersect(varanus, rownames(pca$x))
ptLab <- c(-2, 3.2)
segments(x0 = pca$x[varanus, 1], x1 = rep(ptLab[1], length(varanus)), y0 = pca$x[varanus, 2], y1 = rep(ptLab[2], length(varanus)), lend = 1, lwd = 0.5, col = labLineCol)
text(ptLab[1], ptLab[2], 'Varanus', pos = 3, font = 4, col = 'slateblue1')


ind <- which(rownames(pca$x) %in% lizards)
points(pca$x[ind, 1], pca$x[ind, 2], col ='black', bg = ptCols[ind], pch = 21, cex = ptSizes[ind])

ind <- which(rownames(pca$x) %in% snakes)
points(pca$x[ind, 1], pca$x[ind, 2], col ='black', bg = ptCols[ind], pch = 21, cex = ptSizes[ind])


points(x = grconvertX(seq(0.8, 0.88, length.out = 4), from='npc', to = 'user'), y = rep(grconvertY(0.02, from = 'npc', to = 'user'), 4), cex = scales::rescale(seq(min(breadth, na.rm = TRUE), max(breadth, na.rm = TRUE), length.out = 4), to = c(0.5, 2.5)))
text(x=grconvertX(0.84, from='npc', to = 'user'), y = grconvertY(0.025, from = 'npc', to = 'user'), labels = 'dietary breadth', pos = 3)


# Comparison of diet breadth

plot.new()
plot.window(xlim = c(0,2), ylim = range(breadth, na.rm = TRUE))
axis(1, at = c(0.5, 1.5), labels = c('lizards', 'snakes'), lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l", lwd = 0.5)
mtext('diet breadth', side = 2, line = 2.5, cex = 0.8)

width <- 0.25

violin(breadth[lizards], center = 0.5, width = 0.25, col = 'slateblue1', truncate = TRUE, border = NA)
violin(breadth[snakes], center = 1.5, width = 0.25, col = 'coral2', truncate = TRUE, border = NA)




