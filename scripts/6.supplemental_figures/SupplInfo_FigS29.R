setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

require(ape)
require(sf)
require(rmapshaper)
require(geos)
require(rnaturalearth)

source('./scripts/6.supplemental_figures/violin.R')

alldat <- read.csv('./data/alldat.csv')

tr <- read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
tr <- ladderize(tr, right = TRUE)

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
tr <- rotate(tr, nn)

tr <- read.tree(text = write.tree(tr)) # otherwise, node rotation doesn't change structure

snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']



world <- ne_download(scale = 'small', type = 'land', category = 'physical')
world <- st_geometry(st_as_sf(world))
world <- st_transform(world, crs = '+proj=eqearth')

regionDir <- '~/Dropbox/squamatePhylo/2019/empirical/customRegions'

regionFiles <- list.files(regionDir, pattern='.gpkg$', full.names=TRUE)

## There is an option 1/2 for Central and South America. We will use option 1
regionFiles <- setdiff(regionFiles, grep('option2', regionFiles, value=T))

regionList <- vector('list', length(regionFiles))
names(regionList) <- gsub('\\.gpkg$', '', basename(regionFiles))
for (i in 1:length(regionFiles)) {
	message(i)
	xx <- st_read(regionFiles[i])
	xx <- st_geometry(xx)
	xx <- sfheaders::sf_remove_holes(xx)
	
	xx <- st_transform(xx, "+proj=eqearth")
	xx <- st_combine(xx)
	xx <- as_geos_geometry(xx)
	xx <- geos_unary_union(xx)

	if (!geos_is_valid(xx)) {
		xx <- geos_make_valid(xx)
	}
	
	geos_is_valid(xx)
	
	regionList[[i]] <- xx
	
}

names(regionList) <- gsub('_option\\d', '', names(regionList))

geogRegions <- c('NorthAmerica', 'Palearctic', 'IndoMalayan_mainland', 'SundaShelf', 'CentralAmerica', 'PapuaNewGuinea', 'SouthAmerica', 'WallaceLyddekerTransition', 'SubSaharanAfrica', 'Madagascar', 'Australia', 'NewZealand')

tolower(geogRegions) %in% tolower(names(regionList))


simplifiedPoly <- vector('list', length(geogRegions))
for (i in 1:length(geogRegions)) {
	
	message('\t', i)
	poly <- regionList[[which(tolower(names(regionList)) == tolower(geogRegions[i]))]]
	simplifiedPoly[[i]] <- st_geometry(ms_simplify(st_as_sf(poly)))
}

# combined
combinedRegions <- st_as_sf(data.table::rbindlist(lapply(simplifiedPoly, st_as_sf)))
combinedRegions <- combinedRegions[1:nrow(combinedRegions),]
combinedRegions <- st_geometry(combinedRegions)

combinedRegionsSimplified <- ms_simplify(combinedRegions)


regionCols <- c(
'#8dd3c7',
'#fdbf6f',
'#bebada',
'#fb8072',
'#80b1d3',
'#fdb462',
'#b3de69',
'#fccde5',
'#d9d9d9',
'#bc80bd',
'#ccebc5',
'#ffed6f')


regionCentroids <- st_coordinates(st_centroid(combinedRegionsSimplified))
rownames(regionCentroids) <- geogRegions

regionLabels <- c(	
	NorthAmerica = 'North America',
	Palearctic = 'Palearctic',
	IndoMalayan_mainland = 'Indomalayan mainland',
	SundaShelf = 'Sunda shelf',
	CentralAmerica = 'Central America',
	PapuaNewGuinea = 'New Guinea',
	SouthAmerica = 'South America',
	WallaceLyddekerTransition = 'Wallace Lyddeker transition',
	SubSaharanAfrica = 'Sub-Saharan Africa',
	Madagascar = 'Madagascar',
	Australia = 'Australia',
	NewZealand = 'New Zealand')







mat <- matrix(c(	14, 2, 3, 4, 5,14,
					 6, 1, 1, 1, 1, 7,
					 8, 1, 1, 1, 1, 9,
					14,10,11,12,13,14 
), nrow = 4, ncol = 6, byrow = TRUE)



pdf('~/Downloads/figS29.pdf', width = 8, height = 5)

layout(mat)

par(mar = c(0,0,0,0), oma = c(1,1,1,1))

plot(world, col = gray(0.9), lwd = 0.3, border = gray(0.8))

plot(combinedRegionsSimplified, add = TRUE, lwd = 0.3, col = regionCols)


regionCentroidsDev <- regionCentroids[geogRegions, ]
regionCentroidsDev[,1] <- grconvertX(regionCentroids[geogRegions, 1], from = 'user', to = 'ndc')
regionCentroidsDev[,2] <- grconvertY(regionCentroids[geogRegions, 2], from = 'user', to = 'ndc')

par(mar = c(1,2,1,1))

for (i in 1:length(geogRegions)) {
	
	sp <- alldat[alldat[, paste0('geog_', geogRegions[i])] == 1, 'treename']
	sp <- intersect(sp, tr$tip.label)
	snakeSp <- intersect(sp, snakes)
	lizardSp <- intersect(sp, lizards)
	
	# exclude sea snakes
	seasnakes <- grep('^Hydrophis_|^Laticauda_|^Parahydrophis_|^Hydrelaps_|^Ephalophis_|^Emydocephalus_|^Aipysurus_', snakeSp, value = TRUE)
	snakeSp <- setdiff(snakeSp, seasnakes)
	
	snakeRates <- log(setNames(alldat[, 'meanImputedCLADS'], alldat$treename)[snakeSp])
	lizardRates <- log(setNames(alldat[, 'meanImputedCLADS'], alldat$treename)[lizardSp])
	
	plot.new()
	plot.window(xlim = c(0.5, 2), ylim = range(log(alldat$meanImputedCLADS), na.rm = TRUE))
	#axis(1, lwd = 0, lwd.ticks = 1)
	axis(2, lwd = 0, lwd.ticks = 1, las = 1, cex.axis = 0.65, mgp = c(3, 0.65, 0))
	box(which = "plot", bty = "o")
	mtext('log CLaDS speciation rate', side = 2, line = 1.5, cex = 0.5)
	
	if (length(snakeSp) > 0) {
		violin(snakeRates, center = 1, truncate = TRUE, width = 0.15, col = adjustcolor('coral2', alpha.f = 0.75), border = NA, xpd = NA)
	}
	if (length(lizardSp) > 0) {
		violin(lizardRates, center = 1.5, truncate = TRUE, width = 0.15, col = adjustcolor('slateblue1', alpha.f = 0.75), border = NA, xpd = NA)
	}
	
	# top panels: (0.5, 0)
	
	if (i %in% 1:4) {
		startCoord <- c(0.5, 0)
		titleSide <- 3
	} else if (i %in% 9:12) {
		startCoord <- c(0.5, 1)
		titleSide <- 1
	} else if (i %in% c(5,7)) {
		startCoord <- c(1, 0.5)
		titleSide <- 3
	} else if (i %in% c(6,8)) {
		startCoord <- c(0, 0.5)
		titleSide <- 3
	}
	
	mtext(regionLabels[i], side = titleSide, cex = 0.5)
	
	segments(	x0 = grconvertX(startCoord[1], from = 'npc', to = 'user'), y0 = grconvertY(startCoord[2], from = 'npc', to = 'user'), 
				x1 = grconvertX(regionCentroidsDev[i, 1], from  = 'ndc', to = 'user'),
				y1 = grconvertY(regionCentroidsDev[i, 2], from  = 'ndc', to = 'user'),
				xpd = NA)
	points(grconvertX(regionCentroidsDev[i, 1], from  = 'ndc', to = 'user'), grconvertY(regionCentroidsDev[i, 2], from  = 'ndc', to = 'user'), xpd = NA, pch = 20)
	
}

legend(x = grconvertX(1.2, from = 'npc', 'user'), y = grconvertY(0.7, from = 'npc', 'user'), legend = c('snakes', 'lizards'), fill = adjustcolor(c('coral2', 'slateblue1', alpha.f = 0.75)), bty = 'n', border = NA, xpd = NA, cex = 1)

dev.off()













