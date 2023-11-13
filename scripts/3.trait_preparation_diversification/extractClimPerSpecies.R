# Process CHELSA bioclimatic variables + GMTED2010 elevation (TRI) and record mean values per species.


require(sf)
require(terra)
require(geos)
require(rmapshaper)

# file locations
geogfile <- '~/Dropbox/squamatePhylo/2019/empirical/GARD1.1_dissolved_ranges/modeled_reptiles.shp'
elevFile <- '~/Dropbox/elev_GMTED2010.tif'
masterTaxFile <- '~/Dropbox/Oz_Crown_Ages/taxon-attributes/squamatesCoreTaxonomy_repDB.csv'
bioclimDir <- '~/Dropbox/squamatePhylo/chelsa'
setwd('~/Dropbox/squamatePhylo/2019/empirical')

## varanus
geogfile <- '~/squam/modeled_reptiles.shp'
elevFile <- '~/squam/elev_GMTED2010.tif'
masterTaxFile <- '~/squam/squamatesCoreTaxonomy_repDB.csv'
bioclimDir <- '~/squam/chelsa'
setwd('~/squam')

# -----------------------------

spList <- readRDS('squamGeogPolygons.rds')
spList <- lapply(spList, st_geometry)

# Read in bioclim data

bioclim <- rast(list.files(bioclimDir, pattern = 'bio.+\\.tif$', full.names = TRUE))
names(bioclim) <- gsub('(CHELSA_)(bio\\d+)(_1981-2010_V.2.1)', '\\2', names(bioclim))

### Calculate climatic moisture index 
## Inputs needed are annual precipitation (bio12) and annual PET, both of which we can get from CHELSA.

# annualPrecip <- bioclim[['bio12']]
# petFiles <- list.files('~/chelsa', pattern = '\\.tif', full.names = TRUE)
# petFiles <- petFiles[grep('penman', basename(petFiles))]
# monthlyPET <- rast(petFiles)
# PET <- sum(monthlyPET)

# ind <- annualPrecip < PET
# cmi <- rast(annualPrecip)
# cmi[ind == 1] <- (annualPrecip[ind == 1]/PET[ind == 1]) - 1
# cmi[ind == 0] <- 1 - (PET[ind == 0]/annualPrecip[ind == 0])
# names(cmi) <- "cmi"
# writeRaster(cmi, '~/chelsa/cmi.tif')
cmi <- rast(list.files(bioclimDir, pattern = 'cmi.tif$', full.names = TRUE))

npp <- rast(list.files(bioclimDir, pattern = 'npp', full.names = TRUE))
names(npp) <- 'npp'

# Read in elevation data and calculate terrain ruggedness index (mean of the absolute differences between the value of a cell and the value of its 8 surrounding cells)

elev <- rast(elevFile)
tri <- terrain(elev, v = 'TRI')

# combine into one spatRaster
bioclim <- c(bioclim, cmi, npp, elev, tri)
names(bioclim)[22:23] <- c('elev', 'tri')
bioclim <- bioclim[[c(paste0('bio', 1:19), 'cmi', 'npp', 'elev', 'tri')]]

gc()
# tmpFiles(current=T, orphan=T,old=T, remove=T)

# mask all with landmass
library(rnaturalearth)
land <- ne_download(scale = 50, type = 'land', category = 'physical')
land <- st_as_sf(land)
land <- st_geometry(land)

landras <- rasterize(vect(land), rast(bioclim[[1]]))

seqCol <- colorRampPalette(c('#fef0d9', '#fdcc8a', '#fc8d59', '#e34a33', '#b30000'))
divCol <- colorRampPalette(c('#d7191c', '#fdae61', '#ffffbf', '#abdda4', '#2b83ba'))

npp <- mask(npp, vect(land))
cmi <- mask(cmi, vect(land))
png('npp.png', width = 10, height = 7, units = 'in', res = 400)
plot(npp, col = seqCol(100))
dev.off()

png('cmi.png', width = 10, height = 7, units = 'in', res = 400)
plot(cmi, col = divCol(100), range = c(-1, 1))
dev.off()

# for each species, get the mean value for each bioclim variable across the range
bioclimMat <- matrix(nrow = length(spList), ncol = nlyr(bioclim))
colnames(bioclimMat) <- names(bioclim)
rownames(bioclimMat) <- names(spList)

pb = txtProgressBar(min = 0, max = length(spList), style = 3)
for (i in 1:length(spList)) {
	setTxtProgressBar(pb, i)
# 	if (!all(st_is_valid(spList[[i]]))) stop()
	xx <- extract(bioclim, vect(spList[[i]]), fun = mean, na.rm = TRUE)
	bioclimMat[i, ] <- as.numeric(xx[1, colnames(bioclimMat)])
}
close(pb)

saveRDS(bioclimMat, 'bioclimMat.rds')

for (i in 1:length(spList)) {
	if (!all(st_is_valid(spList[[i]]))) {
		message('\trepairing ', i)
		spList[[i]] <- st_make_valid(spList[[i]])
	}
}


# add on centroid, and min and max latitude and longitude
bioclimMat <- cbind(centroidLong = NA, centroidLat = NA, minLong = NA, maxLong = NA, minLat = NA, maxLat = NA, bioclimMat)
pb = txtProgressBar(min = 0, max = length(spList), style = 3)
for (i in 1:length(spList)) {
	setTxtProgressBar(pb, i)

	combined <- st_combine(spList[[i]])
	if (all(st_is_valid(combined))) {
		xy <- st_coordinates(st_centroid(combined))
	} else {
		xy <- st_coordinates(st_centroid(st_combine(st_simplify(spList[[i]], preserveTopology = TRUE, dTolerance = 0.1))))
	}
	bb <- st_bbox(combined)
	
	bioclimMat[i, 'centroidLong'] <- xy[1]
	bioclimMat[i, 'centroidLat'] <- xy[2]
	
	bioclimMat[i, 'minLong'] <- bb$xmin
	bioclimMat[i, 'maxLong'] <- bb$xmax
	bioclimMat[i, 'minLat'] <- bb$ymin
	bioclimMat[i, 'maxLat'] <- bb$ymax

}
close(pb)

# Convert names to tree names (with imputed)
head(bioclimMat)

# add in geographic range size
rangeSizes <- sapply(spList, function(x) sum(expanse(vect(x), unit = 'km')))
names(rangeSizes) <- names(spList)
bioclimMat <- cbind(bioclimMat, rangeSize = rangeSizes[rownames(bioclimMat)])

head(bioclimMat)

write.csv(bioclimMat, 'climateMat.csv', row.names = TRUE)





# Extract climate data for climate space figure

proj <- '+proj=eqearth'

spList <- readRDS('squamGeogPolygons.rds')
spList <- lapply(spList, st_geometry)


# Read in bioclim data
bioclim <- rast(list.files(bioclimDir, pattern = 'bio.+\\.tif$', full.names = TRUE))
names(bioclim) <- gsub('(CHELSA_)(bio\\d+)(_1981-2010_V.2.1)', '\\2', names(bioclim))

cmi <- rast(list.files(bioclimDir, pattern = 'cmi.tif$', full.names = TRUE))

npp <- rast(list.files(bioclimDir, pattern = 'npp', full.names = TRUE))
names(npp) <- 'npp'

# Read in elevation data and calculate terrain ruggedness index (mean of the absolute differences between the value of a cell and the value of its 8 surrounding cells)

# combine into one spatRaster
bioclim <- c(bioclim, cmi, npp)

# exclude cells that have NA for any layer
qq <- sum(bioclim)


pts <- spatSample(qq, size = 1e4, method = 'random', na.rm = TRUE, as.points = TRUE)

# pts <- st_transform(pts, crs = 4326)

# get sampling of climate
env <- extract(bioclim, pts, xy = TRUE)

pts2 <- st_as_sf(pts)

# get climate information for each species
spEnv <- vector('list', length(spList))
names(spEnv) <- names(spList)
for (i in 1:length(spList)) {	
	if (i %% 100 == 0) message('\t', i)
	spEnv[[i]] <- unlist(st_intersects(spList[[i]], pts2))
	# spEnv[[i]] <- extract(bioclim, vect(st_intersection(pts, spList[[i]])))
}

saveRDS(list(env = env, sp = spEnv), 'climateSpace.rds')







