# Impute SVL from total length
## Pull together multiple datasets containing SVL or total length and tail length, perform some taxonomic reconciliation, and impute SVL based on body/tail length ratios.

require(readxl)
require(ape)
require(nlme)
require(scales)
require(Rphylopars)
require(stringr)

setwd('~/Dropbox/Oz_Crown_Ages')

master <- read.csv('./taxon-attributes/squamatesCoreTaxonomy_repDB.csv')

lengthMass <- read.csv('./trait-data/SVLimputation/lengthAndMass.csv')
ummz <- as.data.frame(read_excel('./trait-data/ummz-morphology/SnakeMaster_Cleaned_addedbodysizes_HW_20160527.xlsx'))
feldman2013 <- as.data.frame(read_excel('./trait-data/SVLimputation/Feldman2013_snake_species_length_and_mass_data.xlsx'))
sheehy2016 <- read.csv('./trait-data/SVLimputation/Sheehy2016/SheehyDat.csv')
lawing2012 <- read.csv('./trait-data/SVLimputation/Lawing2012/Lawing2012_relTailLength.csv')
king1989 <- read.csv('./trait-data/SVLimputation/King1989/king1989.csv')
vertnetfile <- './trait-data/SVLimputation/vertnet_squam_withLength-f59c07e50ae740368267e2c352a540dd.txt'

tree <- read.tree('./final-trees/best_ultrametric_fulltree_ddBD_revision.tre')
tree <- read.tree(text = write.tree(ladderize(tree)))

repDBdir <- './reptileDB/'

# Read in Reptile Database data
acceptedNamesTable <- read.csv(paste0(repDBdir, 'acceptedTaxonTable.csv'), stringsAsFactors = FALSE)
subspeciesTable <- read.csv(paste0(repDBdir, 'subspeciesTable.csv'), stringsAsFactors = FALSE)
synonymTable <- read.csv(paste0(repDBdir, 'synonymTable.csv'), stringsAsFactors = FALSE)

# limit synonyms to post-1980
synonymTable <- synonymTable[which(synonymTable$year >= 1980 | is.na(synonymTable$year)),]

source('~/Dropbox/Oz_Crown_Ages/reptileDB/matchToReptileDB.R')

repDB_acceptedBinomial <- apply(acceptedNamesTable, 1, function(x) paste0(x[1], ' ', x[2]))


# name reconciliation for UMMZ

ummz$Tree_Fullname <- stringr::str_to_title(ummz$Tree_Fullname)

missing <- c()
multi <- c()

for (i in 1:nrow(ummz)) {
	
	tax <- ummz$Tree_Fullname[i]
	
	if (tax %in% master$repDB) {
			
		# we can treat this as a direct match
		ind <- which(master$repDB == tax)
		ummz[i, 'repDB'] <- master[ind, 'repDB']
		
	} else if (tax %in% master$treename) {
		# we can treat this as a direct match
		ind <- which(master$treename == tax)
		if (length(ind) > 1) {
			multi <- c(multi, tax)
		} else {
			ummz[i, 'repDB'] <- master[ind, 'repDB']
		}		
	} else {
		
		if (any(grepl(tax, master$geog))) {
			ind <- grep(tax, master$geog)
			if (any(grepl(',', grep(tax, master$geog, value = TRUE))) & length(ind) > 1) {
				ind <- which(!grepl(',', grep(tax, master$geog, value = TRUE)))
				ind <- grep(tax, master$geog)[ind]
			}
			if (length(ind) > 1) stop()

			ummz[i, 'repDB'] <- master[ind, 'repDB']			

		} else if (tax %in% master$eco) {

			ind <- which(master$eco == tax)
			if (length(ind) > 1) stop()
	
				ummz[i, 'repDB'] <- master[ind, 'repDB']
	
		} else {
			missing <- c(missing, tax)
		}	
	}
}


table(is.na(ummz$repDB))
nomatch <- which(is.na(ummz$repDB))

multi
missing

for (i in 1:nrow(ummz)) {
	
	if (is.na(ummz[i, 'repDB'])) {
		
		tax <- ummz$Tree_Fullname[i]
		match <- matchToReptileDB(gsub('_', ' ', tax), considerSynonyms = TRUE, resolveToSubSpecies = FALSE, returnMultipleMatches = TRUE)
		
		if (is.vector(match)) {
			matched <- gsub('\\s+', '_', match['repdbTaxon'])
			if (matched %in% master$repDB) {
				ummz[i, 'repDB'] <- matched
			} else {
				stop()
			}
		} else {
			# if multiple matches but one is just a genus update, pick that one
			ind <- grep(paste0('\\s', gsub('(.+)_(.+)', '\\2', tax), '$'), match[, 1])
			if (length(ind) == 1) {
				ummz[i, 'repDB'] <- gsub('\\s+', '_', match[ind, 1])
			} else {
				multi <- c(multi, tax)
			}
		}
	}
}

for (i in 1:nrow(ummz)) {
	if (!is.na(ummz[i, 'repDB'])) {
		if (ummz[i, 'repDB'] %in% master$repDB) {
			ummz[i, 'treename'] <- master[which(master$repDB == ummz[i, 'repDB']), 'treename']
		} else if (ummz[i, 'repDB'] %in% master$treename) {
			ummz[i, 'treename'] <- ummz[i, 'repDB']
		} else {
			stop()
		}
	}
}

ummz[which(is.na(ummz$treename)), 'treename'] <- ummz[which(is.na(ummz$treename)), 'repDB']


table(is.na(ummz[, 'repDB']))
ummz[is.na(ummz$repDB),]

for (i in 1:length(unique(ummz$repDB))) {
	
	ind <- which(ummz$repDB == unique(ummz$repDB)[i])
	if (length(unique(ummz[ind, 'Tree_Fullname'])) > 1) stop()
	
}




## Name reconciliation for Feldman and Meiri 2013

feldman2013$Species <- gsub('\\s+', '_', feldman2013$Species)

for (i in 1:nrow(feldman2013)) {
	
	tax <- feldman2013$Species[i]
	
	if (tax %in% master$repDB) {
			
		# we can treat this as a direct match
		ind <- which(master$repDB == tax)
		feldman2013[i, 'repDB'] <- master[ind, 'repDB']
		
	} else if (tax %in% master$treename) {
		# we can treat this as a direct match
		ind <- which(master$treename == tax)
		if (length(ind) > 1) {
			multi <- c(multi, tax)
		} else {
			feldman2013[i, 'repDB'] <- master[ind, 'repDB']
		}		
	} else {
		
		if (any(grepl(tax, master$geog))) {
			ind <- grep(tax, master$geog)
			if (any(grepl(',', grep(tax, master$geog, value = TRUE))) & length(ind) > 1) {
				ind <- which(!grepl(',', grep(tax, master$geog, value = TRUE)))
				ind <- grep(tax, master$geog)[ind]
			}
			if (length(ind) > 1) stop()

			feldman2013[i, 'repDB'] <- master[ind, 'repDB']			

		} else if (tax %in% master$eco) {

			ind <- which(master$eco == tax)
			if (length(ind) > 1) stop()
	
				feldman2013[i, 'repDB'] <- master[ind, 'repDB']
	
		} else {
			missing <- c(missing, tax)
		}	
	}
}

table(is.na(feldman2013$repDB))
nomatch <- which(is.na(feldman2013$repDB))


for (i in 1:nrow(feldman2013)) {
	
	if (is.na(feldman2013[i, 'repDB'])) {
		
		tax <- feldman2013$Species[i]
		match <- matchToReptileDB(gsub('_', ' ', tax), considerSynonyms = TRUE, resolveToSubSpecies = FALSE, returnMultipleMatches = TRUE)
		
		if (is.vector(match)) {
			matched <- gsub('\\s+', '_', match['repdbTaxon'])
			if (matched %in% master$repDB) {
				feldman2013[i, 'repDB'] <- matched
			} else {
				stop()
			}
		} else {
			# if multiple matches but one is just a genus update, pick that one
			ind <- grep(paste0('\\s', gsub('(.+)_(.+)', '\\2', tax), '$'), match[, 1])
			if (length(ind) == 1) {
				feldman2013[i, 'repDB'] <- gsub('\\s+', '_', match[ind, 1])
			} else {
				multi <- c(multi, tax)
			}
		}
	}
}

table(is.na(feldman2013[, 'repDB']))
feldman2013[is.na(feldman2013$repDB), 'Species']

feldman2013[which(feldman2013$Species == 'Plagiopholis_tyani'), 'repDB'] <- 'Plagiopholis_styani'
feldman2013[which(feldman2013$Species == 'Pseudoboa_coronate'), 'repDB'] <- 'Pseudoboa_coronata'
feldman2013[which(feldman2013$Species == 'Chitulia_ornata'), 'repDB'] <- 'Hydrophis_ornatus'
feldman2013[which(feldman2013$Species == 'psammophis_trinasalis'), 'repDB'] <- 'Psammophis_trinasalis'
feldman2013[which(feldman2013$Species == 'Broghammerus_reticulates'), 'repDB'] <- 'Malayopython_reticulatus'

toDrop <- c()
for (i in 1:length(unique(feldman2013$repDB))) {
	
	ind <- which(feldman2013$repDB == unique(feldman2013$repDB)[i])
	if (length(unique(feldman2013[ind, 'Species'])) > 1) {
		toDrop <- c(toDrop, setdiff(feldman2013[ind, 'Species'], feldman2013[ind[1], 'repDB']))
	}
}

feldman2013[which(feldman2013$Species %in% toDrop), 'repDB'] <- NA


for (i in 1:nrow(feldman2013)) {
	if (!is.na(feldman2013[i, 'repDB'])) {
		if (feldman2013[i, 'repDB'] %in% master$repDB) {
			feldman2013[i, 'treename'] <- master[which(master$repDB == feldman2013[i, 'repDB']), 'treename']
		} else if (feldman2013[i, 'repDB'] %in% master$treename) {
			feldman2013[i, 'treename'] <- feldman2013[i, 'repDB']
		} else {
			stop()
		}
	}
}

feldman2013[which(is.na(feldman2013$treename)), 'treename'] <- feldman2013[which(is.na(feldman2013$treename)), 'repDB']


for (i in 1:length(unique(feldman2013$repDB))) {
	ind <- which(feldman2013$repDB == unique(feldman2013$repDB)[i])
	if (length(unique(feldman2013[ind, 'Species'])) > 1) stop()
}


############################################

# Assemble training dataset

## Feldman 2013
### SVL and TL are not necessarily from the same individuals. So we will only keep SVL and TL measurements for species where the mass and source fields are identical, as these are most likely the same individual. 

doubleSp <- names(which((table(feldman2013$treename) == 2) == TRUE))

trainingdat1 <- as.data.frame(matrix(nrow = length(doubleSp), ncol = 6))
colnames(trainingdat1) <- c('sp', 'source', 'sex', 'SVL', 'totalLength', 'ratio')
for (i in 1:length(doubleSp)) {
	svlInd <- which(feldman2013$treename == doubleSp[i] & feldman2013$Measurment == 'SVL')
	tlInd <- which(feldman2013$treename == doubleSp[i] & feldman2013$Measurment == 'TL')
	if (length(svlInd) == 0 | length(tlInd) == 0) stop()
	
	# are mass and source fields the same?
	if (abs(1 - feldman2013[svlInd, 'Mass'] / feldman2013[tlInd, 'Mass']) < 0.01) {
		if (feldman2013[svlInd, 'Source'] == feldman2013[tlInd, 'Source']) {
			#if (abs(1 - feldman2013[svlInd, 'Mass'] / feldman2013[tlInd, 'Mass']) > 0) check <- c(check, svlInd, tlInd)
			trainingdat1[i, 1] <- doubleSp[i]
			trainingdat1[i, 2] <- 'Feldman2013'
			trainingdat1[i, 4] <- feldman2013[svlInd, 'Length']
			trainingdat1[i, 5] <- feldman2013[tlInd, 'Length']	
		}
	}
}

trainingdat1 <- trainingdat1[!is.na(trainingdat1$sp), ]
trainingdat1$ratio <- trainingdat1$SVL / trainingdat1$totalLength

## UMMZ dataset
### Measurements come from the same individuals. 

doubleSp <- ummz[which(!is.na(ummz$SVL) & !is.na(ummz$Tail_length)), 'treename']

table(ummz$Gender)
ummz[which(ummz$Gender == 'F'), 'Gender'] <- 'female'
ummz[which(ummz$Gender == 'M'), 'Gender'] <- 'male'
ummz[which(ummz$Gender == 'Male'), 'Gender'] <- 'male'
ummz[grep('\\?', ummz$Gender), 'Gender'] <- NA
ummz$Gender <- tolower(ummz$Gender)

# in remarks, there are mentions of tail tips missing, etc. We should ignore those entries.
badTails <- grep('tail_cut|tail_?broken|significant damage|too damaged|tail may be broken|tail tip missing|end of tail missing|tail is damaged|tail tip ablated|tail appears to be broken|tail cut|tail tip broken|tail tip may be missing|end of tail may be broken off| tail_might_be_broken|tail may have broken', ummz$Remarks)
ummz <- ummz[setdiff(1:nrow(ummz), badTails), ]

trainingdat2 <- as.data.frame(matrix(nrow = nrow(ummz), ncol = 6))
colnames(trainingdat2) <- c('sp', 'source', 'sex', 'SVL', 'totalLength', 'ratio')
for (i in 1:nrow(ummz)) {
	if (!is.na(ummz[i, 'SVL']) & !is.na(ummz[i, 'Tail_length'])) {
		trainingdat2[i, 1] <- ummz[i, 'treename']
		trainingdat2[i, 2] <- 'UMMZ'
		trainingdat2[i, 3] <- ummz[i, 'Gender']
		trainingdat2[i, 4] <- ummz[i, 'SVL']
		trainingdat2[i, 5] <- ummz[i, 'SVL'] + ummz[i, 'Tail_length']
	}
}

trainingdat2 <- trainingdat2[!is.na(trainingdat2$sp), ]
trainingdat2$ratio <- trainingdat2$SVL / trainingdat2$totalLength

table(trainingdat2$sex)


## Sheehy et al. 2016 dataset
### Contains relative tail length (tail length / total body length) and total body length in cm
### all females

trainingdat3 <- as.data.frame(matrix(nrow = nrow(sheehy2016), ncol = 6))
colnames(trainingdat3) <- c('sp', 'source', 'sex', 'SVL', 'totalLength', 'ratio')
for (i in 1:nrow(sheehy2016)) {
	trainingdat3[i, 1] <- sheehy2016[i, 'treename']
	trainingdat3[i, 2] <- 'Sheehy2016'
	trainingdat3[i, 3] <- 'female'
	trainingdat3[i, 4] <- (sheehy2016[i, 'totalBodyLength_cm'] - sheehy2016[i, 'totalBodyLength_cm'] * sheehy2016[i, 'relativeTailLength']) * 10
	trainingdat3[i, 5] <- sheehy2016[i, 'totalBodyLength_cm'] * 10
}
trainingdat3$ratio <- trainingdat3$SVL / trainingdat3$totalLength

## Lawing et al. 2012
### Contains relative tail length but not total body length or SVL or tail length.
### We will convert relative tail length to SVL/totalLength ratio as 1 - relative tail length.
### gender is not specified, is a mix

trainingdat4 <- as.data.frame(matrix(nrow = nrow(lawing2012), ncol = 6))
colnames(trainingdat4) <- c('sp', 'source', 'sex', 'SVL', 'totalLength', 'ratio')
for (i in 1:nrow(lawing2012)) {
	trainingdat4[i, 1] <- lawing2012[i, 'treename']
	trainingdat4[i, 2] <- 'Lawing2012'
	trainingdat4[i, 3] <- NA
	trainingdat4[i, 6] <- 1 - lawing2012[i, 'relativeTailLength']
}

# some corrections - based on looking up primary literature
## Antillotyphlops_hypomethes should be 196/201
## Pliocercus_elapoides should be 0.62
## Leptophis_depressirostris should be 0.63

trainingdat4[trainingdat4$sp == 'Antillotyphlops_hypomethes', 'ratio'] <- 196/201
trainingdat4[trainingdat4$sp == 'Pliocercus_elapoides', 'ratio'] <- 0.62
trainingdat4[trainingdat4$sp == 'Leptophis_depressirostris', 'ratio'] <- 0.63

## King 1989
### contains measurements for both males and females, for all listed species

trainingdat5 <- as.data.frame(matrix(nrow = nrow(king1989) * 2, ncol = 6))
colnames(trainingdat5) <- c('sp', 'source', 'sex', 'SVL', 'totalLength', 'ratio')
counter <- 1
for (i in 1:nrow(king1989)) {
	trainingdat5[counter, 1] <- king1989[i, 'treename']
	trainingdat5[counter, 2] <- 'King1989'
	trainingdat5[counter, 3] <- 'female'
	trainingdat5[counter, 4] <- king1989[i, 'FSVL']
	trainingdat5[counter, 5] <- king1989[i, 'FSVL'] + king1989[i, 'FTL']
	counter <- counter + 1
	trainingdat5[counter, 1] <- king1989[i, 'treename']
	trainingdat5[counter, 2] <- 'King1989'
	trainingdat5[counter, 3] <- 'male'
	trainingdat5[counter, 4] <- king1989[i, 'MSVL']
	trainingdat5[counter, 5] <- king1989[i, 'MSVL'] + king1989[i, 'MTL']
	counter <- counter + 1
}
trainingdat5$ratio <- trainingdat5$SVL / trainingdat5$totalLength




# ##################################################
# # search vertnet for SVL and tail lengths
## Downloaded vertnet data for squamata and has length = TRUE
## Since there would be a very large number of records to investigate, and because quality control may be lower in this dataset than in the others, we will only consider records for species we are otherwise missing. 


vertnet <- data.table::fread(vertnetfile, data.table = FALSE)
vertnet <- vertnet[vertnet$specificepithet != '', ]
vertnet <- vertnet[- grep('\\.', vertnet$specificepithet), ]
vertnet <- vertnet[- grep('\\sx\\s', vertnet$specificepithet), ]
vertnet <- vertnet[- grep('\\sgroup', vertnet$specificepithet), ]
vertnet <- vertnet[- grep('\\,', vertnet$specificepithet), ]
vertnet$gensp <- paste0(vertnet$genus, '_', vertnet$specificepithet)

# we are actually only interested in snakes
cladeDefs <- read.csv('./trait-data/cladeMemberships.csv')
snakefams <- unique(cladeDefs[cladeDefs$Serpentes == 1, 'family'])
setdiff(vertnet$family, snakefams)
snakefams <- c(snakefams, 'Atractaspididae', 'Dipsadidae', 'Natricidae', 'Hydrophiidae', 'Anilidae', 'Pareatidae')

table(vertnet$family %in% snakefams)
vertnet <- vertnet[vertnet$family %in% snakefams, ]


table(vertnet$gensp %in% master$repDB)

for (i in 1:nrow(vertnet)) {
	
	tax <- vertnet$gensp[i]
	
	if (tax %in% master$repDB) {
			
		# we can treat this as a direct match
		ind <- which(master$repDB == tax)
		vertnet[i, 'repDB'] <- master[ind, 'repDB']
		
	} else if (tax %in% master$treename) {
		# we can treat this as a direct match
		ind <- which(master$treename == tax)
		if (length(ind) > 1) {
			multi <- c(multi, tax)
		} else {
			vertnet[i, 'repDB'] <- master[ind, 'repDB']
		}		
	} else {
		
		if (any(grepl(tax, master$geog))) {
			ind <- grep(tax, master$geog)
			if (any(grepl(',', grep(tax, master$geog, value = TRUE))) & length(ind) > 1) {
				ind <- which(!grepl(',', grep(tax, master$geog, value = TRUE)))
				ind <- grep(tax, master$geog)[ind]
			}
			if (length(ind) > 1) stop()

			vertnet[i, 'repDB'] <- master[ind, 'repDB']			

		} else if (tax %in% master$eco) {

			ind <- which(master$eco == tax)
			if (length(ind) > 1) stop()
	
				vertnet[i, 'repDB'] <- master[ind, 'repDB']
	
		} else {
			missing <- c(missing, tax)
		}	
	}
}

table(is.na(vertnet$repDB))
nomatch <- which(is.na(vertnet$repDB))


for (i in 1:nrow(vertnet)) {
	
	if (is.na(vertnet[i, 'repDB'])) {
		
		tax <- vertnet$gensp[i]
		match <- matchToReptileDB(gsub('_', ' ', tax), considerSynonyms = TRUE, resolveToSubSpecies = FALSE, returnMultipleMatches = TRUE)
		
		if (is.vector(match)) {
			matched <- gsub('\\s+', '_', match['repdbTaxon'])
			if (matched %in% master$repDB) {
				vertnet[i, 'repDB'] <- matched
			} else {
				stop()
			}
		} else {
			# if multiple matches but one is just a genus update, pick that one
			ind <- grep(paste0('\\s', gsub('(.+)_(.+)', '\\2', tax), '$'), match[, 1])
			if (length(ind) == 1) {
				vertnet[i, 'repDB'] <- gsub('\\s+', '_', match[ind, 1])
			} else {
				multi <- c(multi, tax)
			}
		}
	}
}

table(is.na(vertnet$repDB))
nomatch <- which(is.na(vertnet$repDB))

unique(vertnet[is.na(vertnet$repDB), 'gensp'])

vertnet[which(vertnet$gensp == 'Elaphe_obsoleta'), 'repDB'] <- 'Pantherophis_obsoletus'
vertnet[which(vertnet$gensp == 'Charina_orcutti'), 'repDB'] <- 'Lichanura_trivirgata'
vertnet[which(vertnet$gensp == 'Thamnodynastes_strigilis'), 'repDB'] <- 'Thamnodynastes_pallidus'
vertnet[which(vertnet$gensp == 'Imantodes_lentiferas'), 'repDB'] <- 'Imantodes_lentiferus'
vertnet[which(vertnet$gensp == 'Pituophis_melanolecus'), 'repDB'] <- 'Pituophis_melanoleucus'
vertnet[which(vertnet$gensp == 'Lampropeltis_getulus'), 'repDB'] <- 'Lampropeltis_getula'
vertnet[which(vertnet$gensp == 'Crotalus_scutellatus'), 'repDB'] <- 'Crotalus_scutulatus'
vertnet[which(vertnet$gensp == 'Bothriechis_schlegali'), 'repDB'] <- 'Bothriechis_schlegelii'
vertnet[which(vertnet$gensp == 'Tropidonophus_mairii'), 'repDB'] <- 'Tropidonophis_mairii'
vertnet[which(vertnet$gensp == 'Leptodeira_latifasciata'), 'repDB'] <- 'Pseudoleptodeira_latifasciata'
vertnet[which(vertnet$gensp == 'Sibynomorphus_catesbyi'), 'repDB'] <- 'Dipsas_catesbyi'
vertnet[which(vertnet$gensp == 'Tantilla_annulata'), 'repDB'] <- 'Tantilla_supracincta'
vertnet[which(vertnet$gensp == 'Tantilla_gracile'), 'repDB'] <- 'Tantilla_gracilis'
vertnet[which(vertnet$gensp == 'Eryx_miliarius'), 'repDB'] <- 'Eryx_miliaris'
vertnet[which(vertnet$gensp == 'Xenodon_colubrinus'), 'repDB'] <- 'Xenodon_rabdocephalus'
vertnet[which(vertnet$gensp == 'Phrynonax_sulphureus'), 'repDB'] <- 'Spilotes_sulphureus'
vertnet[which(vertnet$gensp == 'Rhinophis_sanguinensis'), 'repDB'] <- 'Rhinophis_sanguineus'
vertnet[which(vertnet$gensp == 'Python_cebae'), 'repDB'] <- 'Python_sebae'
vertnet[which(vertnet$gensp == 'Leptotyphlops_brevicauda'), 'repDB'] <- 'Tricheilostoma_bicolor'
vertnet[which(vertnet$gensp == 'Asthenodipsas_lasgalensis'), 'repDB'] <- 'Asthenodipsas_lasgalenensis'
vertnet[which(vertnet$gensp == 'Phimophis_guianense'), 'repDB'] <- 'Phimophis_guianensis'
vertnet[which(vertnet$gensp == 'Typhlops_diardi'), 'repDB'] <- 'Argyrophis_diardii'

toDrop <- c()
for (i in 1:length(unique(vertnet$repDB))) {
	
	ind <- which(vertnet$repDB == unique(vertnet$repDB)[i])
	if (length(unique(vertnet[ind, 'gensp'])) > 1) {
		toDrop <- c(toDrop, setdiff(vertnet[ind, 'gensp'], vertnet[ind[1], 'repDB']))
	}
}

vertnet[which(vertnet$gensp %in% toDrop), 'repDB'] <- NA


for (i in 1:nrow(vertnet)) {
	if (!is.na(vertnet[i, 'repDB'])) {
		if (vertnet[i, 'repDB'] %in% master$repDB) {
			vertnet[i, 'treename'] <- master[which(master$repDB == vertnet[i, 'repDB']), 'treename']
		} else if (vertnet[i, 'repDB'] %in% master$treename) {
			vertnet[i, 'treename'] <- vertnet[i, 'repDB']
		} else {
			stop()
		}
	}
}

vertnet[which(is.na(vertnet$treename)), 'treename'] <- vertnet[which(is.na(vertnet$treename)), 'repDB']

# limit to those taxa we don't already have covered.
newsp <- setdiff(vertnet$treename, c(trainingdat1$sp, trainingdat2$sp, trainingdat3$sp, trainingdat4$sp, trainingdat5$sp))
newsp <- intersect(newsp, tree$tip.label)

vertnet <- vertnet[vertnet$treename %in% newsp,]


tailInd <- Reduce(union, list(
	grep('tail', vertnet$identificationremarks, ignore.case = TRUE),
	grep('tail', vertnet$occurrenceremarks, ignore.case = TRUE),
	grep('tail', vertnet$dynamicproperties, ignore.case = TRUE)))

svlInd <- Reduce(union, list(
	grep('svl|snout', vertnet$identificationremarks, ignore.case = TRUE),
	grep('svl|snout', vertnet$occurrenceremarks, ignore.case = TRUE),
	grep('svl|snout', vertnet$dynamicproperties, ignore.case = TRUE)))

lengthInd <- Reduce(union, list(
	grep('length', vertnet$identificationremarks, ignore.case = TRUE),
	grep('length', vertnet$occurrenceremarks, ignore.case = TRUE),
	grep('length', vertnet$dynamicproperties, ignore.case = TRUE)))
	
ind <- sort(unique(c(tailInd, svlInd, lengthInd)))


# searching for tail is the only one that was useful

vert2 <- vertnet[ind, c('references', 'family', 'treename', 'gensp', 'lengthinmm', 'lengthtype', 'identificationremarks', 'occurrenceremarks', 'dynamicproperties')]
vert2$SVL <- NA
vert2$tailLength <- NA
vert2$totalLength <- NA
head(vert2)

svlPattern <- 'snout(-|\\s)vent(-|\\s)length|SVL|SV|snoutventlength'
totalLengthPattern <- 'total\\ssize|total|\\bTL\\b|TOT'
tailLengthPattern <- 'tail'

# a couple of fixes. Easier to fix here than to come up with regex
vert2$dynamicproperties <- gsub('mmTL', 'mm TL', vert2$dynamicproperties)
vert2$dynamicproperties <- gsub('l200', '1200', vert2$dynamicproperties)
	
# search for indications of broken or stub tail
badTail <- c(
			grep('stub|broken|tip|regenerated|tail\\sincomplete', vert2[, 'occurrenceremarks'], ignore.case = TRUE),
			grep('stub|broken|tip|regenerated|tail\\sincomplete', vert2[, 'identificationremarks'], ignore.case = TRUE),
			grep('stub|broken|tip|regenerated|tail\\sincomplete', vert2[, 'dynamicproperties'], ignore.case = TRUE))

multicatch <- c()
for (i in 1:nrow(vert2)) {
	
	zz <- vert2[i, 'dynamicproperties']
	zz <- gsub('\\{|\\}|\\"', '', zz)
	zz <- strsplit(zz, ',|; ')[[1]]
	
	svl <- grep(svlPattern, ignore.case = TRUE, zz, value = TRUE)
	totalLength <- grep(totalLengthPattern, ignore.case = TRUE, zz, value = TRUE)
	tailLength <- grep(tailLengthPattern, ignore.case = TRUE, zz, value = TRUE)
	
	svl <- grep('\\d+', svl, value = TRUE)
	totalLength <- grep('\\d+', totalLength, value = TRUE)
	tailLength <- grep('\\d+', tailLength, value = TRUE)
	
	if (length(svl) + length(totalLength) + length(tailLength) == 0) {
		
		zz <- vert2[i, 'occurrenceremarks']
		zz <- gsub('\\{|\\}|\\"', '', zz)
		zz <- strsplit(zz, ',|; ')[[1]]
		if (length(zz) > 1) {
			svl <- grep(svlPattern, ignore.case = TRUE, zz, value = TRUE)
			totalLength <- grep(totalLengthPattern, ignore.case = TRUE, zz, value = TRUE)
			tailLength <- grep(tailLengthPattern, ignore.case = TRUE, zz, value = TRUE)
		} else {
			zz <- str_extract_all(zz, '\\w+\\s?=\\s?\\d+(.\\d+)?')[[1]]
			svl <- grep(svlPattern, ignore.case = TRUE, zz, value = TRUE)
			totalLength <- grep(totalLengthPattern, ignore.case = TRUE, zz, value = TRUE)
			tailLength <- grep(tailLengthPattern, ignore.case = TRUE, zz, value = TRUE)
		}
		
		svl <- grep('\\d+', svl, value = TRUE)
		totalLength <- grep('\\d+', totalLength, value = TRUE)
		tailLength <- grep('\\d+', tailLength, value = TRUE)		
	}
	
	if (length(svl) > 1) stop()
	if (length(totalLength) > 1) stop()
	if (length(tailLength) > 1) stop()
	
	if (length(unlist(str_extract_all(svl, '\\d+(\\.\\d+)?'))) > 1) {
		
		if (length(unlist(str_extract_all(svl, '\\d+(\\.\\d+)?(\\s?(mm|cm)\\s?)?SVL'))) == 1) {
			svl <- str_extract(svl, '\\d+(\\.\\d+)?(\\s?(mm|cm)\\s?)?SVL')
			totalLength <- str_extract(totalLength, '\\d+(\\.\\d+)?(\\s?(mm|cm)\\s?)?total')
		} else if (length(unlist(str_extract_all(svl, 'SVL\\s?\\d+(\\.\\d+)?(\\s?(mm|cm)\\s?)'))) == 1) {
			svl <- str_extract(svl, 'SVL\\s?\\d+(\\.\\d+)?(\\s?(mm|cm)\\s?)')
			totalLength <- str_extract(totalLength, 'Total\\s?\\d+(\\.\\d+)?(\\s?(mm|cm)\\s?)')
		}
	}
	
	if (length(unlist(str_extract_all(totalLength, '\\d+(\\.\\d+)?'))) > 1) {
		totalLength <- str_extract_all(tolower(totalLength), '\\d+(\\.\\d+)?\\s(mm|cm)')[[1]]
	}
	if (length(unlist(str_extract_all(tailLength, '\\d+(\\.\\d+)?'))) > 1) {
		tailLength <- str_extract_all(tolower(tailLength), '\\d+(\\.\\d+)?\\s(mm|cm)')[[1]]
	}
	
	if (length(totalLength) == 1) {
		if (grepl('feet', totalLength)) {
			totalLength <- character(0)
		}
	}
	
	if (length(unlist(str_extract_all(svl, '\\d+(\\.\\d+)?'))) > 1) {
		multicatch <- c(multicatch, i)
	}
	if (length(unlist(str_extract_all(tailLength, '\\d+(\\.\\d+)?'))) > 1) {
		multicatch <- c(multicatch, i)
	}
	if (length(unlist(str_extract_all(totalLength, '\\d+(\\.\\d+)?'))) > 1) {
		multicatch <- c(multicatch, i)
	}
	
	if (!i %in% multicatch) {
		if (length(svl) == 1) {
			vert2[i, 'SVL'] <- as.numeric(str_extract(svl, '\\d+(\\.\\d+)?'))
			if (grepl('cm', svl)) {
				vert2[i, 'SVL'] <- vert2[i, 'SVL'] * 10
			}
		}
	
		if (length(totalLength) == 1) {
			vert2[i, 'totalLength'] <- as.numeric(str_extract(totalLength, '\\d+(\\.\\d+)?'))
			if (grepl('cm', totalLength)) {
				vert2[i, 'totalLength'] <- vert2[i, 'totalLength'] * 10
			}
		}
	
		if (length(tailLength) == 1) {
			vert2[i, 'tailLength'] <- as.numeric(str_extract(tailLength, '\\d+(\\.\\d+)?'))
			if (grepl('cm', tailLength)) {
				vert2[i, 'tailLength'] <- vert2[i, 'tailLength'] * 10
			}
		}
		
		if (i %in% badTail) {
			vert2[i, c('SVL', 'totalLength', 'tailLength')] <- NA
		}
	}
		
	vert2[i,]	
	
}

# a few records were difficult to sort out with regex. 
vert2[multicatch,]
vert2[which(vert2$references == 'http://portal.vertnet.org/o/cm/herps?id=urn-catalog-cm-herps-r334'), c('totalLength', 'tailLength')] <- c(334, 4)

vert2[which(vert2$references == 'http://mczbase.mcz.harvard.edu/guid/MCZ:Herp:R-180315'), c('SVL', 'totalLength')] <- c(546, 233)

vert2[which(vert2$references == 'http://portal.vertnet.org/o/ummz/herps?id=urn-catalog-ummz-herps-175761'), c('SVL', 'tailLength')] <- c(283, 28)

vert2[which(vert2$references == 'http://portal.vertnet.org/o/ummz/herps?id=urn-catalog-ummz-herps-176159'), c('SVL', 'tailLength')] <- c(286, 23)

vert2[which(vert2$references == 'http://portal.vertnet.org/o/ummz/herps?id=urn-catalog-ummz-herps-175763'), c('SVL', 'tailLength')] <- c(255, 27)

# Philothamnus_carinatus: significantly shorter than others in this genus, so we will preemptively drop it. 
vert2 <- vert2[vert2$treename != 'Philothamnus_carinatus',]


vert2[which(vert2$SVL == 0 | vert2$totalLength == 0 | vert2$tailLength == 0), c('SVL', 'totalLength', 'tailLength')] <- NA

# which did not result in measurements?
vert2[which(is.na(vert2$SVL) & is.na(vert2$totalLength) & is.na(vert2$tailLength)),]

# check: if tail length greater or equal to SVL? Probably should have been total length
# check: is total length less than SVL? Probably should have been tail length.
vert2[which(vert2[, 'tailLength'] >= vert2[, 'SVL']),]

vert2[which(vert2[, 'totalLength'] <= vert2[, 'SVL']),]

# We will assume that these are in error and that tail and total length were switched.
zz <- which(vert2[, 'tailLength'] >= vert2[, 'SVL'])
vert2[zz, 'totalLength'] <- vert2[zz, 'tailLength']
vert2[zz, 'tailLength'] <- NA

zz <- which(vert2[, 'totalLength'] <= vert2[, 'SVL'])
vert2[zz, 'tailLength'] <- vert2[zz, 'totalLength']
vert2[zz, 'totalLength'] <- NA


vert2 <- vert2[ - which(is.na(vert2$SVL) & is.na(vert2$totalLength) & is.na(vert2$tailLength)),]

# we need 2 out of 3
NAcount <- apply(vert2[, c('SVL', 'totalLength', 'tailLength')], 1, function(x) sum(is.na(x)))
vert2[which(NAcount == 2),]
vert2 <- vert2[NAcount < 2,]

trainingdat6 <- as.data.frame(matrix(nrow = nrow(vert2), ncol = 6))
colnames(trainingdat6) <- c('sp', 'source', 'sex', 'SVL', 'totalLength', 'ratio')
for (i in 1:nrow(vert2)) {
	trainingdat6[i, 1] <- vert2[i, 'treename']
	trainingdat6[i, 2] <- 'VERTNET'
	trainingdat6[i, 3] <- NA
	if (is.na(vert2[i, 'SVL'])) {
		if (is.na(vert2[i, 'totalLength']) & is.na(vert2[i, 'tailLength'])) stop()
		trainingdat6[i, 4] <- vert2[i, 'totalLength'] - vert2[i, 'tailLength']
	} else {
		trainingdat6[i, 4] <- vert2[i, 'SVL']
	}
	if (is.na(vert2[i, 'totalLength'])) {
		if (is.na(vert2[i, 'SVL']) | is.na(vert2[i, 'tailLength'])) stop()
		trainingdat6[i, 5] <- vert2[i, 'SVL'] + vert2[i, 'tailLength']
	} else {
		trainingdat6[i, 5] <- vert2[i, 'totalLength']
	}
}

trainingdat6$ratio <- trainingdat6$SVL / trainingdat6$totalLength

## Combine datasets

# what taxa are unique to each dataset?
setdiff(trainingdat1$sp, c(trainingdat2$sp, trainingdat3$sp, trainingdat4$sp, trainingdat5$sp))
setdiff(trainingdat2$sp, c(trainingdat1$sp, trainingdat3$sp, trainingdat4$sp, trainingdat5$sp))
setdiff(trainingdat3$sp, c(trainingdat2$sp, trainingdat1$sp, trainingdat4$sp, trainingdat5$sp))
setdiff(trainingdat4$sp, c(trainingdat2$sp, trainingdat3$sp, trainingdat1$sp, trainingdat5$sp))
setdiff(trainingdat5$sp, c(trainingdat2$sp, trainingdat3$sp, trainingdat1$sp, trainingdat4$sp))
setdiff(trainingdat6$sp, c(trainingdat1$sp, trainingdat2$sp, trainingdat3$sp, trainingdat4$sp, trainingdat5$sp))

# As vertnet is the only non-curated dataset, we will only incorporate for species we are otherwise missing.

trainingdat6 <- trainingdat6[trainingdat6$sp %in% setdiff(trainingdat6$sp, c(trainingdat1$sp, trainingdat2$sp, trainingdat3$sp, trainingdat4$sp, trainingdat5$sp)), ]

alldat <- rbind(trainingdat1, trainingdat2, trainingdat3, trainingdat4, trainingdat5, trainingdat6)

alldat <- alldat[order(alldat$sp), ]

# any where total length is the same or less than SVL?
## would indicate some problem or data entry error
which(alldat$totalLength <= alldat$SVL)

head(alldat)
table(alldat$source)
sort(table(alldat$sp[which(duplicated(alldat$sp))]))

# if we were to restrict ourselves to males or females, how many species would we have?
length(unique(alldat[which(alldat$sex == 'male'), 'sp']))
length(unique(alldat[which(alldat$sex == 'female'), 'sp']))

## about 300

setdiff(unique(alldat[which(alldat$sex == 'male'), 'sp']), unique(alldat[which(alldat$sex == 'female'), 'sp']))
setdiff(unique(alldat[which(alldat$sex == 'female'), 'sp']), unique(alldat[which(alldat$sex == 'male'), 'sp']))

# without restricting by sex, we have 1247 species.

# For those species that have multiple individuals, we will keep the larger one. 
## If there are males and females, we will keep the largest female. 
## If there are males or females and some not specified, we will ignore the unspecified sex.

# how much variation is there in SVL/TL ratio within species?
ratioRange <- lapply(split(alldat, alldat$sp), function(x) range(x$ratio))
ratioRange[order(sapply(ratioRange, diff))]
table(sapply(ratioRange, diff) > 0.1)

# remove clearly problematic records
alldat[alldat$ratio <= 0.4,]
alldat <- alldat[alldat$ratio > 0.4,]



# Preferentially keep values from published datasets. 
# order of preference:
# Sheehy
# King, Lawing, Feldman
# UMMZ
# VERTNET 
# If we have multiple records for a species, we will pick the one closest to the median of the SVL/TL ratio.

keepIndTraining <- c()
for (i in 1:length(unique(alldat$sp))) {
	
	taxon <- unique(alldat$sp)[i]
	ind <- which(alldat$sp == taxon)
	ind2 <- which(alldat$sp == taxon & alldat$source == 'Sheehy2016')
	ind3 <- which(alldat$sp == taxon & alldat$source %in% c('Feldman2013', 'King1989', 'Lawing2012'))
	ind4 <- which(alldat$sp == taxon & alldat$source == 'UMMZ')
	ind5 <- which(alldat$sp == taxon & alldat$source == 'VERTNET')
	if (length(ind) > 1) {
		# if (anyNA(alldat[ind, 'SVL'])) stop()
		
		# Sheehy 2016
		if (length(ind2) > 0) {
			
			medianVal <- median(alldat[ind2, 'ratio'])
			hasSVL <- intersect(ind2, which(alldat$sp == taxon & !is.na(alldat$SVL)))
			
			# for training dataset
			if (length(hasSVL) > 0) {
				keepIndTraining <- c(keepIndTraining, hasSVL[which.min(abs(medianVal - alldat[hasSVL, 'ratio']))])
			} else {
				keepIndTraining <- c(keepIndTraining, ind2[which.min(abs(medianVal - alldat[ind2, 'ratio']))])
			}

			
		# King 1989, Feldman 2013, Lawing 2012	
		} else if (length(ind3) > 0) {

			medianVal <- median(alldat[ind3, 'ratio'])
			hasSVL <- intersect(ind3, which(alldat$sp == taxon & !is.na(alldat$SVL)))
			
			if (length(hasSVL) > 0) {
				keepIndTraining <- c(keepIndTraining, hasSVL[which.min(abs(medianVal - alldat[hasSVL, 'ratio']))])
			} else {
				keepIndTraining <- c(keepIndTraining, ind3[which.min(abs(medianVal - alldat[ind3, 'ratio']))])
			}
			
			
		# UMMZ
		} else if (length(ind4) > 0) {
			
			medianVal <- median(alldat[ind4, 'ratio'])
			hasSVL <- intersect(ind4, which(alldat$sp == taxon & !is.na(alldat$SVL)))
			
			# for training dataset
			if (length(hasSVL) > 0) {
				keepIndTraining <- c(keepIndTraining, hasSVL[which.min(abs(medianVal - alldat[hasSVL, 'ratio']))])
			} else {
				keepIndTraining <- c(keepIndTraining, ind4[which.min(abs(medianVal - alldat[ind4, 'ratio']))])
			}
			


		# VERTNET
		} else if (length(ind5) > 0) {
			
			medianVal <- median(alldat[ind5, 'ratio'])
			hasSVL <- intersect(ind5, which(alldat$sp == taxon & !is.na(alldat$SVL)))
			
			# for training dataset
			if (length(hasSVL) > 0) {
				keepIndTraining <- c(keepIndTraining, hasSVL[which.min(abs(medianVal - alldat[hasSVL, 'ratio']))])
			} else {
				keepIndTraining <- c(keepIndTraining, ind5[which.min(abs(medianVal - alldat[ind5, 'ratio']))])
			}
			

						
		}		

	} else if (length(ind) == 1) {
		
		# for training dataset
		keepIndTraining <- c(keepIndTraining, ind)
		
	}
}


trainingdat <- alldat[keepIndTraining, ]
testingdat <- trainingdat[!is.na(trainingdat$SVL),]

nrow(trainingdat)
nrow(testingdat)

anyDuplicated(trainingdat$sp)
anyDuplicated(testingdat$sp)

trainingdat[trainingdat$sp %in% trainingdat$sp[duplicated(trainingdat$sp)],]

all(keepIndTesting %in% keepIndTraining)
setdiff(keepIndTesting, keepIndTraining)

# several taxa are not in tree. Drop them.
setdiff(trainingdat[, 'sp'], tree$tip.label)
trainingdat <- trainingdat[trainingdat$sp %in% tree$tip.label,]

setdiff(testingdat[, 'sp'], tree$tip.label)
testingdat <- testingdat[testingdat$sp %in% tree$tip.label,]

# what is the phylogenetic distribution of the taxa in the training set?
snakeTree <- extract.clade(tree, getMRCA(tree, c('Namibiana_labialis', 'Coluber_constrictor')))
snakeTree <- read.tree(text = write.tree(ladderize(snakeTree)))


# ----------------------------------------------------------------------
# Impute SVL/TL ratio and then estimate SVL from TL using those ratios


# quick test to understand how phylopars works. 
df <- trainingdat[, c('sp', 'ratio')]
colnames(df) <- c('species', 'ratio')

df <- cbind.data.frame(species = snakeTree$tip.label, hasDat = snakeTree$tip.label %in% trainingdat$sp, ratio = NA)
df$ratio <- setNames(trainingdat$ratio, trainingdat$sp)[df$species]

xx <- phylopars(trait_data = df, tree = snakeTree, model = 'BM')
head(xx$anc_recon)


# do leave-one-out cross-validation, imputing SVL/TL ratio on all snakes at once


imputedRatio <- rep(NA, nrow(trainingdat))
names(imputedRatio) <- trainingdat$sp

for (i in 1:nrow(trainingdat)) {
	
	message('\t', i)
		
	focalSp <- trainingdat[i, 'sp']
		
	df <- cbind.data.frame(species = snakeTree$tip.label, hasDat = snakeTree$tip.label %in% trainingdat$sp, ratio = NA)
	df$ratio <- setNames(trainingdat$ratio, trainingdat$sp)[df$species]
	
	df[df$species == focalSp, 'ratio'] <- NA
	
	xx <- phylopars(trait_data = df, tree = snakeTree, model = 'BM')
	imputedRatio[i] <- xx$anc_recon[focalSp, 'ratio']
}



		
# instead of phylogenetic imputation, just take a local mean SVL/TL ratio (this is just for comparison)
imputedRatio2 <- rep(NA, nrow(trainingdat))
for (i in 1:nrow(trainingdat)) {
		
	spVec <- byNode(tree, trainingdat[i, 'sp'], nTaxa = 11, taxonSet = unique(trainingdat$sp), criterion = 'taxa')
	spVec <- setdiff(spVec, trainingdat[i, 'sp'])
	
	imputedRatio2[i] <- mean(trainingdat[match(spVec, trainingdat$sp), 'ratio'])
}	


identical(names(imputedRatio), trainingdat$sp)
anyNA(imputedRatio)

pdf('./trait-data/SVLimputation/phyloImputationRatio.pdf', width = 12, height = 7)
par(mfrow=c(1,2))

plot(trainingdat$ratio, imputedRatio, xlab = 'true ratio', ylab = 'imputed ratio', asp = 1)
abline(a = 0, b = 1, col = 'red')
title(main = paste0('phylogenetic imputation \nof ratio, R2 = ', round(summary(lm(imputedRatio ~ trainingdat$ratio))$adj.r.squared, 3)))

summary(lm(imputedRatio ~ trainingdat$ratio))
summary(lm(imputedRatio2 ~ trainingdat$ratio))

# now with imputed ratios, we can get SVL from TL. 
plot(testingdat$SVL, testingdat$totalLength * imputedRatio[testingdat$sp], xlab = 'true SVL', ylab = 'imputed SVL', asp = 1)
abline(a = 0, b = 1, col = 'red')
title(main = paste0('SVL calculated via imputation \nof SVL/TL ratio, R2 = ', round(summary(lm(testingdat$totalLength * imputedRatio[testingdat$sp] ~ testingdat$SVL))$adj.r.squared, 3)))

summary(lm(testingdat$totalLength * imputedRatio[testingdat$sp] ~ testingdat$SVL))

dev.off()

# Do full imputation using all available data

df <- trainingdat[, c('sp', 'ratio')]
colnames(df) <- c('species', 'ratio')

df <- cbind.data.frame(species = snakeTree$tip.label, hasDat = snakeTree$tip.label %in% trainingdat$sp, ratio = NA)
df$ratio <- setNames(trainingdat$ratio, trainingdat$sp)[df$species]

xx <- phylopars(trait_data = df, tree = snakeTree, model = 'BM')
head(xx$anc_recon)

allRatios <- setNames(xx$anc_recon[snakeTree$tip.label, 'ratio'], snakeTree$tip.label)





# plot true and imputed values, and then plot outlined bars showing imputed values for values that we know. 

pdf(height=55, width=12, file = "./trait-data/SVLimputation/snake_SVL-TL_ratio_imputationCrossValidation.pdf")
layout(matrix(1:2, nrow = 1, ncol = 2), widths = c(90, 10))
qq <- plot.phylo(snakeTree, show.tip.label = TRUE, edge.width = 0.2, cex = 0.15)
xval <- 135

segments(x0 = xval, y0 = 1:Ntip(snakeTree), x1 = 10 + xval, y1 =  1:Ntip(snakeTree), xpd = NA, lend = 1, col = gray(0.95))

cols <- rep('cornflowerblue', Ntip(snakeTree))
cols[which(!snakeTree$tip.label %in% trainingdat$sp)] <- 'coral'
table(cols)

segments(x0 = xval, y0 = 1:Ntip(snakeTree), x1 = rescale(xx$anc_recon[snakeTree$tip.label, 'ratio'], from = c(0, 1), to = c(0, 10)) + xval, y1 =  1:Ntip(snakeTree), xpd = NA, lend = 1, col = cols)
text(xval + 4, Ntip(snakeTree) + 8, labels = 'SVL / TL', cex = 0.5, xpd = NA)
text(xval + 4, Ntip(snakeTree) + 5, labels = 'imputed', cex = 0.3, xpd = NA, col = 'coral')
text(xval + 4, Ntip(snakeTree) + 3, labels = 'known', cex = 0.3, xpd = NA, col = 'cornflowerblue')


imputedInd <- which(snakeTree$tip.label %in% trainingdat$sp)
# match(trainingdat$sp, snakeTree$tip.label)

arrows(x0 = xval, y0 = imputedInd, x1 = rescale(imputedRatio[na.omit(match(snakeTree$tip.label, trainingdat$sp))], from = c(0, 1), to = c(0, 10)) + xval, y1 = imputedInd, xpd = NA, lend = 1, code = 2, angle = 90, lwd = 0.1, length = 0.01)


# legend()

shift <- 20

# error in imputed ratio -- imputed ratio / true ratio
segments(x0 = xval + shift, x1 = xval + shift, y0 = 0, y1 = Ntip(snakeTree), col = gray(0.85), lend = 1, lwd = 0.2)
segments(x0 = xval + shift, y0 = 1:Ntip(snakeTree), x1 = rescale(setNames(imputedRatio/trainingdat$ratio, trainingdat$sp)[snakeTree$tip.label], from = c(0.5, 1.5), to = c(-5, 5)) + xval + shift, y1 =  1:Ntip(snakeTree), xpd = NA, lend = 1, col = 'black')

# axis
segments(x0 = xval + shift - 5, x1 = xval + shift + 5, y0 = -2, y1 = -2, lwd = 0.25, lend = 1)
segments(x0 = xval + shift + c(-5, 0, 5), y0 = -2 - 1, y1 = -2 + 1, lwd = 0.25, lend = 1)
text(xval + shift + c(-5, 0, 5), y = -2.5, labels = c(0.5, 1, 1.5), pos = 1, cex = 0.25)

segments(x0 = xval + shift - 5, x1 = xval + shift + 5, y0 = Ntip(snakeTree) + 1, y1 = Ntip(snakeTree) + 1, lwd = 0.25, lend = 1)
segments(x0 = xval + shift + c(-5, 0, 5), y0 = Ntip(snakeTree) + 1 - 1, y1 = Ntip(snakeTree) + 1 + 1, lwd = 0.25, lend = 1)
text(xval + shift + c(-5, 0, 5), y = Ntip(snakeTree) + 3, labels = c(0.5, 1, 1.5), cex = 0.25)


text(xval + shift, Ntip(snakeTree) + 8, labels = 'imputed / true', cex = 0.5, xpd = NA)

dev.off()

 
# ------------------------------

# Now take the Feldman et al. 2015 dataset of squamate body sizes and estimate SVL for all taxa for which we only have total length, given the SVL/TL ratios imputed above.

head(lengthMass)

lengthMass$completeSVL <- NA

for (i in 1:nrow(lengthMass)) {
	
	if (!is.na(lengthMass[i, 'maxLength'])) {
		
		# if this taxon is not in tree, but there is SVL, then copy it over.
		if (lengthMass[i, 'lengthType'] == 'SVL') {
			lengthMass[i, 'completeSVL'] <- lengthMass[i, 'maxLength'] 
		}
		
		# for taxa with total length and that are not imputed in the tree, estimate SVL
		if (lengthMass[i, 'lengthType'] == 'TL' & lengthMass[i, 'dataType'] != 'imputed') {
			lengthMass[i, 'completeSVL'] <- allRatios[lengthMass[i, 'treename']] * lengthMass[i, 'maxLength']
		}
	}
}


# calculate elongation index by approximating squamates as a cylinder with length as the species' SVL and volume as the species' mass.

# mass (volume) = pi * r^2 * SVL
# r = sqrt(mass / (SVL * pi))
# elongation index = length / diameter = SVL / 2r


lengthMass$elongationIndex = lengthMass$completeSVL / (2*sqrt(lengthMass$mass / (pi * lengthMass$completeSVL)))



pdf('./trait-data/SVLimputation/elongationIndex.pdf', width = 6, height = 100)
par(fig=c(0,0.9,0,1),mar=c(3,0,2,0))
plot(tree, edge.color=8, show.tip.label=TRUE, cex = 0.1, edge.width = 0.25)
nodelabels('', getMRCA(tree, snakeTree$tip.label), frame = 'circle', cex = 0.5)

par(fig=c(0.9,0.97,0,1),new=TRUE, mar=c(3,0,2,0))
plot(setNames(log(lengthMass$completeSVL), lengthMass$treename)[tree$tip.label], 1:ape::Ntip(tree), xaxt='n', yaxt='n', bty='n', type='l', lwd = 0.25)
abline(v=mean(log(lengthMass$completeSVL), na.rm = TRUE), lty=2, col=2, lwd = 0.25)
text(x= mean(log(lengthMass$completeSVL), na.rm = TRUE), y = Ntip(tree), labels = 'log SVL', pos = 3, xpd = NA, cex = 0.25)

par(fig=c(0.97,1,0,1),new=TRUE, mar=c(3,0,2,0))
plot(setNames(lengthMass$elongationIndex, lengthMass$treename)[tree$tip.label], 1:ape::Ntip(tree), xaxt='n', yaxt='n', bty='n', type='l', lwd = 0.25)
text(x= mean(lengthMass$elongationIndex, na.rm = TRUE), y = Ntip(tree), labels = 'elongation\nratio', pos = 3, xpd = NA, cex = 0.25)

dev.off()




write.csv(lengthMass, './trait-data/SVLimputation/lengthMassComplete.csv', row.names = FALSE)
