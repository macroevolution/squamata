# Stacked barplot showing total proportion utilized by snakes vs lizards


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




pdf('~/Downloads/figS18.pdf', width = 14, height = 7)

par(mar = c(1, 4, 2, 1), oma = c(0,0,1,0), mfrow = c(1,2))

# PC 1
nPC <- 1

ind <- order(pca$rotation[, nPC], decreasing = TRUE)

plot.new()
plot.window(xlim = c(1, length(ind)), ylim = max(abs(pca$rotation[, nPC])) * c(-1,1))
axis(2)
mtext('PC 1 loadings', side = 2, line = 2.5)

width <- 0.3

for (i in 1:length(ind)) {
	snakeProp <- dietProportions[intersect(rownames(dietProportions), snakes), ind[i]]
	lizardProp <- dietProportions[intersect(rownames(dietProportions), lizards), ind[i]]
	snakePercent <- sum(snakeProp) / (sum(snakeProp) + sum(lizardProp))
	lizardPercent <- sum(lizardProp) / (sum(snakeProp) + sum(lizardProp))
	dir <- ifelse(pca$rotation[, nPC][ind[i]] > 0, 1, -1)
	rect(xleft = i - width, xright = i + width, ybottom = 0, ytop = rescale(snakePercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])), col = 'coral2', lwd = 0.5)
	rect(xleft = i - width, xright = i + width, ybottom = rescale(snakePercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])), ytop = rescale(snakePercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])) + rescale(lizardPercent, from = c(0,1), to = c(0, pca$rotation[, nPC][ind[i]])), col = 'slateblue1', lwd = 0.5)	
}

pos <- which(pca$rotation[ind, nPC] > 0)
neg <- which(pca$rotation[ind, nPC] < 0)

text(x = pos, y = -0.01, labels = gsub('_', ' ', names(pos)), pos = 2, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)
text(x = neg, y = +0.01, labels = gsub('_', ' ', names(neg)), pos = 4, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)

legend('top', legend = c('snakes', 'lizards'), fill = c('coral2', 'slateblue1'), bty = 'n')



# PC 2
nPC <- 2

ind <- order(pca$rotation[, nPC], decreasing = TRUE)

plot.new()
plot.window(xlim = c(1, length(ind)), ylim = max(abs(pca$rotation[, nPC])) * c(-1,1))
axis(2)
mtext('PC 2 loadings', side = 2, line = 2.5)

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

mtext('Relative cumulative proportions of prey type', side = 3, font = 2, cex = 1.2, outer = TRUE, line = -1)

dev.off()


