# climate space

setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

require(ape)


alldat <- read.csv('./data/alldat.csv')

tr <- read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
tr <- ladderize(tr)

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
tr <- rotate(tr, nn)

tr <- read.tree(text = write.tree(tr)) # otherwise, node rotation doesn't change structure

snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']


env <- readRDS('./data/3.trait_preparation_diversification/climateSpace.rds')
spList <- env[[2]]
env <- env[[1]]
env <- env[complete.cases(env),]

datCol <- 2:22

envpca <- prcomp(env[, datCol], center = TRUE, scale. = TRUE)

# we can project species data into this pca space with predict()
snakeInd <- intersect(snakes, names(spList))
snakeInd <- unique(unlist(spList[snakeInd]))
snakeProj <- predict(envpca, env[env$ID %in% snakeInd, datCol])

lizardInd <- intersect(lizards, names(spList))
lizardInd <- unique(unlist(spList[lizardInd]))
lizardProj <- predict(envpca, env[env$ID %in% lizardInd, datCol])


varLabels <- c(
	bio1 = 'annual mean temp',
	bio2 = 'mean diurnal temp range',
	bio3 = 'isothermality',
	bio4 = 'temp seasonality',
	bio5 = 'max temp warmest month',
	bio6 = 'min temp coldest month',
	bio7 = 'temp annual range',
	bio8 = 'mean temp wettest quarter',
	bio9 = 'mean temp driest quarter',
	bio10 = 'mean temp warmest quarter',
	bio11 = 'mean temp coldest quarter',
	bio12 = 'annual precip',
	bio13 = 'precip wettest month',
	bio14 = 'precip driest month',
	bio15 = 'precip seasonality',
	bio16 = 'precip wettest quarter',
	bio17 = 'precip driest quarter',
	bio18 = 'precip warmest quarter',
	bio19 = 'precip coldest quarter',
	cmi = 'climatic moisture index',
	npp = 'net primary productivity'
)


letterHeight <- -0.5

# Figure comparing occupancy of climate space by lizards and snakes, with PC loadings

pdf('~/Downloads/figS21.pdf', width = 8, height = 8.5)

layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
par(mar = c(2,4,2,2))

plot.new()
plot.window(xlim = range(envpca$x[,1]), ylim = range(envpca$x[,2]))
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")
mtext(paste0('climate PC 1 (', round(summary(envpca)$importance[2,1]*100,1), '%)'), side = 1, line = 2.5, cex = 0.9)
mtext(paste0('climate PC 2 (', round(summary(envpca)$importance[2,2]*100,1), '%)'), side = 2, line = 2.5, cex = 0.9)

points(envpca$x[, 1:2], pch = 20, col = gray(0.8))
points(lizardProj[, 1:2], col = 'slateblue1', lwd = 1.2)

mtext('A', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

# polygon(lizardProj[chull(lizardProj[, 1:2]), 1:2], border = 'slateblue1', lwd = 0.5)
# polygon(snakeProj[chull(snakeProj[, 1:2]), 1:2], border = 'coral2', lwd = 0.5)

plot.new()
plot.window(xlim = range(envpca$x[,1]), ylim = range(envpca$x[,2]))
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")
mtext(paste0('climate PC 1 (', round(summary(envpca)$importance[2,1]*100,1), '%)'), side = 1, line = 2.5, cex = 0.9)
mtext(paste0('climate PC 2 (', round(summary(envpca)$importance[2,2]*100,1), '%)'), side = 2, line = 2.5, cex = 0.9)

points(envpca$x[, 1:2], pch = 20, col = gray(0.8))
points(snakeProj[, 1:2], col = 'coral2', lwd = 1.2)

mtext('B', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

# PC 3/4
plot.new()
plot.window(xlim = range(envpca$x[,3]), ylim = range(envpca$x[,4]))
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")
mtext(paste0('climate PC 3 (', round(summary(envpca)$importance[2,3]*100,1), '%)'), side = 1, line = 2.5, cex = 0.9)
mtext(paste0('climate PC 4 (', round(summary(envpca)$importance[2,4]*100,1), '%)'), side = 2, line = 2.5, cex = 0.9)

points(envpca$x[, 3:4], pch = 20, col = gray(0.8))
points(lizardProj[, 3:4], col = 'slateblue1', lwd = 1.2)

mtext('C', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

legend('topright', legend = c('snakes', 'lizards'), col = c('coral2', 'slateblue1'), pt.lwd = 1.2, bty = 'n', cex = 1, pch = 1)

plot.new()
plot.window(xlim = range(envpca$x[,3]), ylim = range(envpca$x[,4]))
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")
mtext(paste0('climate PC 3 (', round(summary(envpca)$importance[2,3]*100,1), '%)'), side = 1, line = 2.5, cex = 0.9)
mtext(paste0('climate PC 4 (', round(summary(envpca)$importance[2,4]*100,1), '%)'), side = 2, line = 2.5, cex = 0.9)

points(envpca$x[, 3:4], pch = 20, col = gray(0.8))
points(snakeProj[, 3:4], col = 'coral2', lwd = 1.2)

mtext('D', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

#### PC LOADINGS

# PC 1
nPC <- 1

ind <- order(envpca$rotation[, nPC], decreasing = TRUE)

plot.new()
plot.window(xlim = c(1, length(ind)), ylim = max(abs(envpca$rotation[, 1:2])) * c(-1,1))
axis(2, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0))
mtext('PC 1 loadings', side = 2, line = 2.5, cex = 0.9)

width <- 0.3

rect(xleft = (1:length(ind)) - width, xright = (1:length(ind)) + width, ybottom = 0, ytop = envpca$rotation[, nPC][ind], col = gray(0.9), lwd = 0.5)

pos <- which(envpca$rotation[ind, nPC] > 0)
neg <- which(envpca$rotation[ind, nPC] < 0)

text(x = pos, y = -0.01, labels = gsub('_', ' ', varLabels[names(pos)]), pos = 2, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)
text(x = neg, y = +0.01, labels = gsub('_', ' ', varLabels[names(neg)]), pos = 4, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)

mtext('E', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)


# PC 2
nPC <- 2

ind <- order(envpca$rotation[, nPC], decreasing = TRUE)

plot.new()
plot.window(xlim = c(1, length(ind)), ylim = max(abs(envpca$rotation[, 1:2])) * c(-1,1))
axis(2, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0))
mtext('PC 2 loadings', side = 2, line = 2.5, cex = 0.9)

width <- 0.3

rect(xleft = (1:length(ind)) - width, xright = (1:length(ind)) + width, ybottom = 0, ytop = envpca$rotation[, nPC][ind], col = gray(0.9), lwd = 0.5)

pos <- which(envpca$rotation[ind, nPC] > 0)
neg <- which(envpca$rotation[ind, nPC] < 0)

text(x = pos, y = -0.01, labels = gsub('_', ' ', varLabels[names(pos)]), pos = 2, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)
text(x = neg, y = +0.01, labels = gsub('_', ' ', varLabels[names(neg)]), pos = 4, srt = 90, xpd = NA, cex = 0.75, offset = 0, xpd = NA)

mtext('F', adj = -0.1, padj = letterHeight, font = 2, xpd = NA, cex = 1.2)

dev.off()


# which points are unoccupied?
# data(World, package='tmap')
# World <- st_geometry(World)
# plot(World, lwd = 0.5, col = gray(0.9), border = NA)
# points(env[setdiff(env$ID, union(snakeInd, lizardInd)), c('x', 'y')], cex = 0.5)

length(snakeInd) / nrow(env)
length(lizardInd) / nrow(env)

length(setdiff(snakeInd, lizardInd)) / nrow(env)
length(setdiff(lizardInd, snakeInd)) / nrow(env)


