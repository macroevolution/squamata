# Plot counts of biomes per family
## This would show that families are widespread

# https://geospatial.tnc.org/datasets/TNC::terrestrial-ecoregions/about

require(sf)
require(ape)
require(geos)

setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

alldat <- read.csv('./data/alldat.csv')

tr <- read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
tr <- ladderize(tr)

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
tr <- rotate(tr, nn)

tr <- read.tree(text = write.tree(tr)) # otherwise, node rotation doesn't change structure

snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']


setwd('~/Dropbox/squamatePhylo/2019/empirical')

# GARD range polygons 
spList <- readRDS('squamGeogPolygons.rds')
spList <- lapply(spList, st_geometry)



eco <- st_read('~/Downloads/Terrestrial_Ecoregions.gdb.zip')
eco <- st_transform(eco, st_crs(spList[[1]]))

apply(st_drop_geometry(eco), 2, function(x) length(unique(x)))



table(st_is_valid(eco))

eco <- st_make_valid(eco)

eco_geos <- as_geos_geometry(eco)

spList <- lapply(spList, as_geos_geometry)

for (i in 1:length(spList)) {
	if (!all(geos_is_valid(spList[[i]]))) {
		message('\trepairing ', i)
		spList[[i]] <- geos_make_valid(spList[[i]])
		if (!all(geos_is_valid(spList[[i]]))) {
			message('\t\t', i, ' still invalid.')
		}
	}
}

spInd <- vector('list', length(spList))
names(spInd) <- names(spList)

pb <- txtProgressBar(max = length(spList), style = 3)
for (i in 1:length(spList)) {
	
	setTxtProgressBar(pb, i)
	if (length(spList[[i]]) == 1) {
		spInd[[i]] <- which(geos_intersects(spList[[i]], eco_geos) == TRUE)	
	} else {
		spInd[[i]] <- unique(unlist(lapply(1:length(spList[[i]]), function(x) which(geos_intersects(spList[[i]][x], eco_geos) == TRUE))))
	}
}
close(pb)

# saveRDS(spInd, '~/Downloads/spInd_biomes.rds')

spInd <- readRDS('~/Downloads/spInd_biomes.rds')

fams <- unique(alldat[, 'family'])


# biomes
famMat <- matrix(0, nrow = length(unique(st_drop_geometry(eco)[, 'WWF_MHTNAM'])), ncol = length(fams))
rownames(famMat) <- sort(unique(st_drop_geometry(eco)[, 'WWF_MHTNAM']))
colnames(famMat) <- sort(fams)
for (i in 1:length(fams)) {
	
	sp <- alldat[alldat$family == fams[i], 'treename']
	sp <- intersect(sp, names(spList))
	ind <- sort(unique(unlist(spInd[sp])))
	xx <- unique(st_drop_geometry(eco)[ind, 'WWF_MHTNAM'])
	famMat[xx, fams[i]] <- famMat[xx, fams[i]] + 1
	
}

# remove Rock/Ice
famMat <- famMat[- which(rownames(famMat) == 'Rock and Ice'), ]

rowSums(famMat)
colSums(famMat)

# realms
famMatRealms <- matrix(0, nrow = length(unique(st_drop_geometry(eco)[, 'WWF_REALM2'])), ncol = length(fams))
rownames(famMatRealms) <- sort(unique(st_drop_geometry(eco)[, 'WWF_REALM2']))
colnames(famMatRealms) <- sort(fams)
for (i in 1:length(fams)) {
	
	sp <- alldat[alldat$family == fams[i], 'treename']
	sp <- intersect(sp, names(spList))
	ind <- sort(unique(unlist(spInd[sp])))
	xx <- unique(st_drop_geometry(eco)[ind, 'WWF_REALM2'])
	famMatRealms[xx, fams[i]] <- famMatRealms[xx, fams[i]] + 1
	
}

rowSums(famMatRealms)
colSums(famMatRealms)

# drop Antarctica
famMatRealms <- famMatRealms[- which(rownames(famMatRealms) == 'Antarctic'), ]


famTreeTips <- sapply(fams, function(x) intersect(tr$tip.label, alldat[alldat$family == x, 'treename'])[1])
famTree <- keep.tip(tr, famTreeTips)
famTree$tip.label <- names(famTreeTips)[match(famTree$tip.label, famTreeTips)]
famTree <- ladderize(famTree, right = FALSE)
famTree <- read.tree(text = write.tree(famTree, file = ''))

famMat <- famMat[, famTree$tip.label]
famMat <- famMat[order(rowSums(famMat), decreasing = TRUE),]

famMatRealms <- famMatRealms[, famTree$tip.label]
famMatRealms <- famMatRealms[order(rowSums(famMatRealms), decreasing = TRUE),]




pdf('~/Downloads/figS22.pdf', width = 14, height = 7.5)
par(mar = c(1, 10, 20, 1))
plot.phylo(famTree, direction = 'upwards', cex = 0.75)

refs <- grconvertY(c(0.4, 0.7), from = 'ndc', 'user')
refs <- seq(from = refs[1], to = refs[2], length.out = nrow(famMat))

for (i in 1:ncol(famMat)) {
	
	cols <- rep(gray(0.9), nrow(famMat))
	cols[which(famMat[,i] == 1)] <- 'dark green'
	points(x = rep(i, length(refs)), y = refs, pch = 22, bg = cols, lwd = 0.25, xpd = NA, cex = 2)
}

text(x = rep(0.5, length(refs)), y = refs, gsub(' and ', ' & ', rownames(famMat)), pos = 2, cex = 0.5, xpd = NA)


refs <- grconvertY(c(0.75, 0.88), from = 'ndc', 'user')
refs <- seq(from = refs[1], to = refs[2], length.out = nrow(famMatRealms))

for (i in 1:ncol(famMatRealms)) {
	
	cols <- rep(gray(0.9), nrow(famMatRealms))
	cols[which(famMatRealms[,i] == 1)] <- 'dark green'
	points(x = rep(i, length(refs)), y = refs, pch = 22, bg = cols, lwd = 0.25, xpd = NA, cex = 2)
}

text(x = rep(0.5, length(refs)), y = refs, rownames(famMatRealms), pos = 2, cex = 0.5, xpd = NA)


dev.off()