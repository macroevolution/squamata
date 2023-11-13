setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

library(scales)

alldat <- read.csv('./data/alldat.csv')

snakes <- alldat[which(alldat[, 'clade_Serpentes'] == 1), 'treename']
#snakes <- snakes[snakes %in% tree$tip.label]
lizards <- alldat[which(alldat[, 'clade_Serpentes'] == 0), 'treename']
#lizards <- lizards[lizards %in% tree$tip.label]


# stacked barplot showing total proportion utilized by snakes vs lizards

dietProportions <- read.csv('./data/3.trait_preparation_diversification/dietModeling/diet-proportions-phylo.csv', row.names = 1)
dietProportions <- data.matrix(dietProportions)

dietProportionsStan <- sweep(log(dietProportions), 1, rowMeans(log(dietProportions))) # clr

pca = prcomp(dietProportionsStan)

# diet breadth
breadth <- setNames(alldat$dietBreadth, alldat$treename)
breadth <- breadth[!is.na(breadth)]



lizardInd <- which(rownames(pca$x) %in% lizards)
snakeInd <- which(rownames(pca$x) %in% snakes)

ptCols <- rep('slateblue1', nrow(pca$x))
ptCols[rownames(pca$x) %in% snakes] <- 'coral2'
table(ptCols)

ptSizes <- scales::rescale(breadth[rownames(pca$x)], to = c(0.5, 2.5))

pdf('~/Downloads/figS19.pdf', width = 10, height = 12)
layout(matrix(c(1,1,1,1,2,3), nrow = 3, ncol = 2, byrow = TRUE))
par(mar = c(4,4,8,7))
plot.new()
plot.window(xlim = range(pca$x[,3]), ylim = range(pca$x[,4]))

axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "o")
mtext(paste0('diet PC 3 (', round(summary(pca)$importance[2,3]*100,1), '%)'), side = 1, line = 2.5)
mtext(paste0('diet PC 4 (', round(summary(pca)$importance[2,4]*100,1), '%)'), side = 2, line = 2.5)


points(pca$x[lizardInd, 3], pca$x[lizardInd, 4], col ='black', bg = ptCols[lizardInd], pch = 21, cex = ptSizes[lizardInd])
points(pca$x[snakeInd, 3], pca$x[snakeInd, 4], col ='black', bg = ptCols[snakeInd], pch = 21, cex = ptSizes[snakeInd])

points(x = grconvertX(seq(0.02, 0.1, length.out = 4), from='npc', to = 'user'), y = rep(grconvertY(0.03, from = 'npc', to = 'user'), 4), cex = scales::rescale(seq(min(breadth, na.rm = TRUE), max(breadth, na.rm = TRUE), length.out = 4), to = c(0.5, 2.5)))
text(x=grconvertX(0.06, from='npc', to = 'user'), y = grconvertY(0.04, from = 'npc', to = 'user'), labels = 'dietary breadth', pos = 3)

legend('topleft', legend = c('snakes', 'lizards'), fill = c('coral2', 'slateblue1'), bty = 'n', cex = 1.5)

mtext('lizards and snakes, PC 3 & 4', outer = TRUE, line = -3, side = 3, font = 2, cex = 1)
mtext('A', adj = 0, padj = -1, font = 2, cex = 1.2)

#plot.new()
#plot.window(xlim = range(pca$x[,3]), ylim = range(0, max(density(pca$x[lizardInd, 3])$y)))
dens <- density(pca$x[snakeInd, 3])
dens <- cbind(dens$x, dens$y)
dens[,2] <- rescale(dens[,2], to = c(grconvertY(1, from = 'npc', 'user'), 3.8))
dens[dens[,1] > grconvertX(1, from = 'npc', 'user'), 1] <- grconvertX(1, from = 'npc', 'user')
dens[dens[,1] < grconvertX(0, from = 'npc', 'user'), 1] <- grconvertX(0, from = 'npc', 'user')
polygon(dens, col = adjustcolor('coral2', alpha.f = 0.75), border = NA, xpd = NA)
dens <- density(pca$x[lizardInd, 3])
dens <- cbind(dens$x, dens$y)
dens[,2] <- rescale(dens[,2], to = c(grconvertY(1, from = 'npc', 'user'), 3.8))
dens[dens[,1] > grconvertX(1, from = 'npc', 'user'), 1] <- grconvertX(1, from = 'npc', 'user')
dens[dens[,1] < grconvertX(0, from = 'npc', 'user'), 1] <- grconvertX(0, from = 'npc', 'user')
polygon(dens, col = adjustcolor('slateblue1', alpha.f = 0.75), border = NA, xpd = NA)

# plot.new()
# plot.window(xlim = range(pca$x[,4]), ylim = range(0, max(density(pca$x[lizardInd, 4])$y)))
dens <- density(pca$x[snakeInd, 4])
dens <- cbind(dens$x, dens$y)
dens2 <- dens
dens2[,1] <- dens[,2]
dens2[,2] <- dens[,1]
dens2[,1] <- rescale(dens2[,1], to = c(grconvertX(1, from = 'npc', 'user'), 5))
dens2[dens2[,2] > grconvertY(1, from = 'npc', 'user'), 2] <- grconvertY(1, from = 'npc', 'user')
dens2[dens2[,2] < grconvertY(0, from = 'npc', 'user'), 2] <- grconvertY(0, from = 'npc', 'user')
polygon(dens2, col = adjustcolor('coral2', alpha.f = 0.75), border = NA, xpd = NA)
dens <- density(pca$x[lizardInd, 4])
dens <- cbind(dens$x, dens$y)
dens2 <- dens
dens2[,1] <- dens[,2]
dens2[,2] <- dens[,1]
dens2[,1] <- rescale(dens2[,1], to = c(grconvertX(1, from = 'npc', 'user'), 5))
dens2[dens2[,2] > grconvertY(1, from = 'npc', 'user'), 2] <- grconvertY(1, from = 'npc', 'user')
dens2[dens2[,2] < grconvertY(0, from = 'npc', 'user'), 2] <- grconvertY(0, from = 'npc', 'user')
polygon(dens2, col = adjustcolor('slateblue1', alpha.f = 0.75), border = NA, xpd = NA)


par(mar = c(1,4,1,1))

# PC 3
nPC <- 3

ind <- order(pca$rotation[, nPC], decreasing = TRUE)

plot.new()
plot.window(xlim = c(1, length(ind)), ylim = max(abs(pca$rotation[, nPC])) * c(-1,1))
axis(2)
mtext('PC 3 loadings', side = 2, line = 2.5)

width <- 0.3

for (i in 1:length(ind)) {
	snakeProp <- dietProportions[intersect(rownames(dietProportions), snakes), ind[i]]
	lizardProp <- dietProportions[intersect(rownames(dietProportions), lizards), ind[i]]
	snakePercent <- sum(snakeProp) / (sum(snakeProp) + sum(lizardProp))
	lizardPercent <- sum(lizardProp) / (sum(snakeProp) + sum(lizardProp))
	rect(xleft = i - width, xright = i + width, ybottom = 0, ytop = rescale(snakePercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])), col = 'coral2', lwd = 0.5)
	rect(xleft = i - width, xright = i + width, ybottom = rescale(snakePercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])), ytop = rescale(snakePercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])) + rescale(lizardPercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])), col = 'slateblue1', lwd = 0.5)	
}

pos <- which(pca$rotation[ind, nPC] > 0)
neg <- which(pca$rotation[ind, nPC] < 0)

text(x = pos, y = -0.01, labels = gsub('_', ' ', names(pos)), pos = 2, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)
text(x = neg, y = +0.01, labels = gsub('_', ' ', names(neg)), pos = 4, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)

mtext('B', adj = 0, padj = -1, font = 2, cex = 1.2)


# PC 4
nPC <- 4

ind <- order(pca$rotation[, nPC], decreasing = TRUE)

plot.new()
plot.window(xlim = c(1, length(ind)), ylim = max(abs(pca$rotation[, nPC])) * c(-1,1))
axis(2)
mtext('PC 4 loadings', side = 2, line = 2.5)

width <- 0.3

for (i in 1:length(ind)) {
	snakeProp <- dietProportions[intersect(rownames(dietProportions), snakes), ind[i]]
	lizardProp <- dietProportions[intersect(rownames(dietProportions), lizards), ind[i]]
	snakePercent <- sum(snakeProp) / (sum(snakeProp) + sum(lizardProp))
	lizardPercent <- sum(lizardProp) / (sum(snakeProp) + sum(lizardProp))
	rect(xleft = i - width, xright = i + width, ybottom = 0, ytop = rescale(snakePercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])), col = 'coral2', lwd = 0.5)
	rect(xleft = i - width, xright = i + width, ybottom = rescale(snakePercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])), ytop = rescale(snakePercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])) + rescale(lizardPercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])), col = 'slateblue1', lwd = 0.5)	
}

pos <- which(pca$rotation[ind, nPC] > 0)
neg <- which(pca$rotation[ind, nPC] < 0)

text(x = pos, y = -0.01, labels = gsub('_', ' ', names(pos)), pos = 2, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)
text(x = neg, y = +0.01, labels = gsub('_', ' ', names(neg)), pos = 4, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)

# mtext('Relative cumulative proportions of prey type', side = 3, font = 2, cex = 1.2, outer = TRUE, line = -1)

dev.off()


