


require(pbapply)
require(ape)
require(sf)

setwd('~/Dropbox/squamatePhylo/2019/empirical')

outfile <- '~/Dropbox/squamatePhylo/2019/empirical/masterTranslationTable.csv'

repDBdir <- '~/Dropbox/Oz_Crown_Ages/reptileDB/'

geogFile <- '~/Dropbox/squamatePhylo/2019/empirical/GARD1.1_dissolved_ranges/modeled_reptiles.shp'
traitFile <- '~/Dropbox/Oz_Crown_Ages/trait-data/Meiri2018_GEB/geb12773-sup-0001-appendixs1.csv'

# read in metadata table for phylogeny
taxonTableFile <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/fasta_March2020/metatable2.csv'
if (file.exists(taxonTableFile)) {
	taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
} else {
	taxonTableFile <- '~/squam2020/metatable2.csv'
	taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
}
head(taxonTable)

# Corrections
taxonTable[which(taxonTable$genus %in% c('Myersophis', 'Oxyrhabdium', 'Hologerrhum', 'Cyclocorus' )), 'family'] <- 'Lamprophiidae'
taxonTable[which(taxonTable$genus == 'Xylophis'), 'family'] <- 'Pareidae'
taxonTable[which(taxonTable$genus == 'Kladirostratus'), 'family'] <- 'Lamprophiidae'

taxonTable[which(taxonTable$genus %in% c("Hologerrhum", "Cyclocorus",  "Myersophis", "Oxyrhabdium")), 'family'] <- 'Cyclocoridae'
taxonTable[which(taxonTable$genus == 'Micrelaps'), 'family'] <- 'Micrelapidae'
taxonTable[which(taxonTable$genus == 'Buhoma'), 'family'] <- 'Elapidae'









pyronTaxonomy <- read.csv('~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/species_family.csv')
pyronTaxonomy$Genus <- sapply(strsplit(pyronTaxonomy$Species, '\\s+|_'), function(x) x[[1]])

# read in tree 
contree <- read.tree('~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/tier1/fulltree/bestTrees/fulltree_default_con_1_raxmlOpt.raxml.bestTree')

# quick adjustment
contree$tip.label[grep('gehyra_cysp', contree$tip.label, ignore.case = TRUE)] <- 'Gehyra_CYsp'

Ntip(contree)
length(unique(taxonTable$repdbTaxon))

# get the actual list of taxa in the tree
treeTaxa <- sort(contree$tip.label)
treeTaxa <- setdiff(treeTaxa, 'Sphenodon_punctatus')
treeTaxa <- gsub('_', ' ', treeTaxa)

# Read in Roll et al. metadata for geographic ranges
geog <- st_read(geogFile)
table(st_drop_geometry(geog)$Group)
geog <- geog[which(!geog$Group %in% c('croc', 'turtle')),]
table(st_drop_geometry(geog)$Group)
geogTaxa <- st_drop_geometry(geog)$Binomial
rm(geog)
geogTaxa <- gsub('^\\s+|\\s+$', '', geogTaxa)
geogTaxa <- sort(unique(gsub('\\s+', ' ', geogTaxa)))

# ecological data
eco <- read.csv(traitFile)
table(eco$Family)
ecoTaxa <- sort(unique(gsub('\\s+', ' ', eco$Binomial)))
ecoTaxa <- gsub('^\\s+|\\s+$', '', ecoTaxa)
ecoTaxa <- setdiff(ecoTaxa, '')

# Read in Reptile Database data
acceptedNamesTable <- read.csv(paste0(repDBdir, 'acceptedTaxonTable.csv'), stringsAsFactors = FALSE)
subspeciesTable <- read.csv(paste0(repDBdir, 'subspeciesTable.csv'), stringsAsFactors = FALSE)
synonymTable <- read.csv(paste0(repDBdir, 'synonymTable.csv'), stringsAsFactors = FALSE)

# limit synonyms to post-1980
synonymTable <- synonymTable[which(synonymTable$year >= 1980 | is.na(synonymTable$year)),]

source('~/Dropbox/Oz_Crown_Ages/reptileDB/matchToReptileDB.R')

repDB_acceptedBinomial <- apply(acceptedNamesTable, 1, function(x) paste0(x[1], ' ', x[2]))

# are any of the tree taxa not in the accepted taxa list from Rep DB?
zz <- taxonTable[which(taxonTable$repdbTaxon %in% setdiff(treeTaxa, repDB_acceptedBinomial)), ]
lapply(split(zz, zz$repdbTaxon), function(x) x[1, c('ncbiTaxon', 'ncbiTaxonID', 'repdbTaxon', 'reason')])


#############################################
#############################################

########################################
# ULTIMATELY THESE NAME CHANGES TO TAXA IN THE TREE WILL BE IN TREE, BUT FOR NOW, MAKE THESE CHANGES HERE

# What changes to make to tree taxon names?
# Xenochrophis flavipunctatus renamed to Fowlea flavipunctatus
# Xenochrophis piscator renamed to Fowlea piscator
# Ptychozoon kuhli renamed to Gekko kuhli.
# Chilomeniscus stramineus renamed to Sonora straminea.

treeTaxa[grep('Xenochrophis flavipunctatus', treeTaxa)] <- 'Fowlea flavipunctatus'
treeTaxa[grep('Xenochrophis piscator', treeTaxa)] <- 'Fowlea piscator'
treeTaxa[grep('Ptychozoon kuhli', treeTaxa)] <- 'Gekko kuhli'
treeTaxa[grep('Chilomeniscus stramineus', treeTaxa)] <- 'Sonora straminea'

##########################################

# We will use accepted taxon names from Reptile DB as the master column. 
# We will then associate taxa in the tree, geog and trait data to those master taxa.
# In cases where the taxon is not present in Reptile DB (for instance, where we opted to keep a taxon present in NCBI genbank), we will still include those, but they will have NA as reptile DB.
# In cases, where our tree taxon matches to multiple Reptile DB taxa, we will list them multiple times. 

masterTable <- matrix(nrow = length(repDB_acceptedBinomial), ncol = 6)
colnames(masterTable) <- c('repDB', 'tree', 'geog', 'eco', 'family', 'genus')
masterTable <- as.data.frame(masterTable)
masterTable$repDB <- repDB_acceptedBinomial
addons <- list()

resList <- replicate(length(repDB_acceptedBinomial), NA, simplify = FALSE)
names(resList) <- repDB_acceptedBinomial

singleMatches <- c()
for (i in 1:length(treeTaxa)) {
	
	tax <- treeTaxa[i]
	
	#message('\t', tax)
	
	# is this taxon present in the reptile BD?
	if (tax %in% repDB_acceptedBinomial) {
		
		# we can treat this as a direct match
		ind <- which(masterTable$repDB == tax)
		resList[[ind]] <- c(setdiff(resList[[ind]], NA), tax)
	
	} else {
		
		# taxon is not immediately found in reptile DB
		
		# try matching
		match <- matchToReptileDB(tax, considerSynonyms = TRUE, resolveToSubSpecies = FALSE, returnMultipleMatches = TRUE)
		
		# was a match found via exact matching?
		if ('authority' %in% colnames(match)) {
			
			if (nrow(match) == 1) {
				# single match
				stop('single match accepted')
			}
			
			ind <- which(masterTable$repDB %in% match$taxon)
			for (j in ind) {
				resList[[j]] <- c(setdiff(resList[[j]], NA), tax)
			}
		}
		
		# was a match found via fuzzy matching and/or synonymy?
		if (!'authority' %in% colnames(match)) {				
			if (!all(is.na(match['repdbTaxon']))) {
				# something was found
				if (is.null(dim(match))) {
					# a single match was found
					singleMatches <- c(singleMatches, tax)
					resList[[which(masterTable$repDB %in% match['repdbTaxon'])]] <- tax
					
				} else if (!is.null(dim(match))) {
					# multiple matches found
					
					# in this case, we will fill in the table for multiple Rep DB taxa

					ind <- which(masterTable$repDB %in% match$taxon)
					for (j in ind) {
						resList[[j]] <- c(setdiff(resList[[j]], NA), tax)
					}
				}
			} else {
				# No matches found. We will add these to the table.
				addons
				new <- rep(NA, ncol(masterTable))
				names(new) <- colnames(masterTable)
				new['tree'] <- tax
				addons <- c(addons, list(new))		
			}
		}		
	}
}

which(lengths(resList) > 1)

## Adelphicos quadrivirgatus maps to:
## 		- Adelphicos quadrivirgatus -> Adelphicos quadrivirgatum, Adelphicos newmanorum
##			SOLUTION: A. quadrivirgatus could either be referring to Adelphicos quadrivirgatus newmanorum or to Adelphicos quadrivirgatum. We will leave it such that the tree taxon A. quadrivirgatus is listed for both. It will sort itself out in subsequent analyses. 
## 			NEW SOLUTION: Set Adelphicos_newmanorum to NA for the tree. The tree tip will point to A. quadrivirgatus

## Anatololacerta anatolica and Anatololacerta oertzeni both map as: 
##	 	- Anatololacerta anatolica -> A. anatolica
## 		- Anatololacerta oertzeni -> A. anatolica, A. budaki, A. danfordi, A. pelasgiana
## 		- Anatololacerta budaki -> A. budaki
## 		- Anatololacerta danfordi -> A. danfordi
## 		- Anatololacerta pelasgiana -> A. pelasgiana, A. budaki
## 			SOLUTION: keep A. oertzeni as additional taxon. A. pelasgiana stays as is.

## Anolis rodriguezii and Anolis rodriguezii: 
##	 	- Anolis rodriguezii is the correct name. Need to drop Anolis rodriguezi from tree.

## Anolis whitemani -> Anolis breslini, Anolis saxatilis
## 			Could be Anolis saxatilis (synonymized in 2019), or Anolis whitemani breslini, which was elevated to full species in 1980s. Given that A. saxilitis synonymy is very new, and A. breslini option is quite a bit older, most likely that A. whitemani -> A. saxatilis. 
##				SOLUTION: Map A. whitemani to A. saxatilis. Remove A. whitemani from A. breslini.

## Aspidoscelis sexlineatus and Aspidoscelis sexlineata: 
##	 		SOLUTION: Aspidoscelis sexlineatus is the correct name. Need to drop Aspidoscelis sexlineata from tree.

## Chilomeniscus stramineus
## 		Chilomeniscus stramineus -> Sonora cincta, Sonora fasciata, Sonora straminea
## 			Chilomeniscus stramineus split into C. cincta, fasciata, straminea by Grismer et al. 2002.
## 			C. stramineus was listed as synonym to Sonora cincta in 2002.
## 			A majority of the loci are from Davis-Rabosky et al. 2016, and are labeled as Sonora straminea.
## 				SOLUTION: Can pretty confidently rename to Sonora straminea.

## Cryptoblepharus pannosus and Cryptoblepharus carnabyi: 
##	 	- Cryptoblepharus pannosus -> same, accepted name
##		- Cryptoblepharus carnabyi -> C. australis, C. pannosus, C. plagiocephalus
## 			SOLUTION: C. carnabyi appears to have been split up. ReptileDB: HORNER 2007 also synonymized Cryptoblepharus carnabyi STORR 1976 with C. plagiocephalus.
##				Since we have it and can't confidently remap it, we will leave it in as an additional species.
## 				Drop C. carnabyi from C. plagiocephalus
##				Drop C. carnabyi from C. australis

## Eumeces schneiderii and Eumeces schneideri
## 			SOLUTION: Eumeces schneiderii is the correct name. Need to drop Eumeces schneideri from tree.

## Hemidactylus anamallensis maps to:
## 		- Hemidactylus anamallensis -> Dravidogecko anamallensis, Dravidogecko douglasadamsi, Dravidogecko janakiae, Dravidogecko meghamalaiensis, Dravidogecko septentrionalis, Dravidogecko tholpalli
## 			This appears to be a recent splitting of the taxon. We do not have enough information to assign. 
## 			SOLUTION: Simply list H. anamallensis for all Dravidogecko. Will need to be sorted out in subsequent analyses. But change is so new, that likely will be irrelevant. 
## 			NEW SOLUTION: set Dravidogecko anamallensis, Dravidogecko douglasadamsi, Dravidogecko janakiae, Dravidogecko meghamalaiensis, Dravidogecko septentrionalis, Dravidogecko tholpalli to NA in tree

## Hemidactylus triedrus and Hemidactylus subtriedrus and Hemidactylus whitakeri
## 		- Hemidactylus triedrus -> H. lankae, H. sahgali, H. triedrus, H. whitakeri
## 		- Hemidactylus subtriedrus -> H. kangerensis, H. triedrus, H. whitakeri
##		- Hemidactylus whitakeri -> H. whitakeri
## 			SOLUTION: 	H. triedrus stays as H. triedrus
##						H. subtriedrus mapped to H. triedrus, therefore we should drop H. subtriedrus from tree.
##						H. whitakeri stays as whitakeri
## 				"Hemidactylus subtriedrus has been synonymized with H. triedrus by Mirza et al. 2018, after serveral previous authors had suspected that they are synonymous."

## Hypsiglena ochrorhynchus, Hypsiglena ochrorhyncha, Hypsiglena unaocularus
## 		- Hypsiglena unaocularus -> same, accepted
## 		- Hypsiglena ochrorhynchus -> same, accepted
##		- Hypsiglena ochrorhyncha -> H. ochrorhynchus, H. unaocularus
## 			SOLUTION: Hypsiglena ochrorhyncha cannot be easily resolved. We should just drop Hypsiglena ochrorhyncha

## Lycodon subannulatus and Dryocalamus davisonii
## 		- Lycodon subannulatus -> same, accepted
## 		- Dryocalamus davisonii -> Lycodon davisonii, Lycodon subannulatus
## 			SOLUTION: Map Dryocalamus davisonii to Lycodon davisonii

## Macrovipera schweizeri, Macrovipera lebetina, Macrovipera lebetinus
## 		- Macrovipera schweizeri -> same, accepted
## 		- Macrovipera lebetinus -> same, accepted
## 		- Macrovipera lebetina -> Macrovipera lebetinus, Macrovipera schweizeri
## 			SOLUTION: Macrovipera lebetina should be mapped to Macrovipera lebetinus (but not renamed)

## Madascincus polleni, Madascincus intermedius
## 		- Madascincus polleni -> M. polleni, M. miafina
## 		- Madascincus intermedius -> M. miafina, M. polleni
## 		- Madascincus miafina -> accepted
## 			SOLUTION: M. intermedius should map to M. polleni, but M. polleni already in tree and accounted for. So add as additional taxon.

## Magliophis exiguus and Magliophis exiguum
## 		Magliophis exiguum should be Magliophis exiguus. Should we drop exiguum from tree? Probably.

## Mochlus sundevallii Mochlus sundevalli
## 		Should be Mochlus sundevallii, should drop Mochlus sundevalli from tree.

## Sinonatrix percarinatus and Sinonatrix percarinata
## 		Should be Sinonatrix percarinatus, should drop Sinonatrix percarinata from tree.


## Xenochrophis flavipunctatus and Xenochrophis piscator map as: 
##	 	- Xenochrophis flavipunctatus -> Fowlea flavipunctatus, Fowlea melanzostus
##		- Xenochrophis piscator -> Fowlea piscator, Fowlea melanzostus
## 			SOLUTION: UNCLEAR! Will call X. flavipunctatus Fowlea flavipunctatus, and X. piscator Fowlea piscator
## 						We will just assume that Fowlea melanzostus is not present in our tree.

## Gloydius intermedius and Gloydius saxatilis (not an accepted name) both map as: 
##	 	- Gloydius intermedius -> G. intermedius
##	 	- Gloydius saxatilis -> G. intermedius, G. shedaoensis (also an accepted taxon)
## 			SOLUTION: We will add Gloydius saxitilis as an additional name, since it has its own tip in the tree, and also has geog range info.
##						Alternative of collapsing it with G. intermedius would also be possible. 
## 						From RepDB: Synonymy partly after KHALIKOV & ANANJEVA (pers. comm.). G. saxatilis is not recognized by GUMPRECHT et al. 2004. ORLOV & BARABANOV 1999 synonymized G. saxatilis with G. intermedius.

## Kladirostratus acutus and Psammophylax acutus map as: 
##	 	- Kladirostratus acutus -> Kladirostratus acutus
##		- Psammophylax acutus -> Kladirostratus acutus, Kladirostratus togoensis
##			SOLUTION: Cannot resolve based on names alone. But probably pretty safe to call Psammophylax acutus Kladirostratus acutus
## 						So we will drop Psammophylax acutus from the phylogeny.

## Ptychozoon lionotum
## 		- Ptychozoon lionotum -> Gekko lionotum, Gekko cicakterbang
## 			Gekko cicakterbang is a new taxon (2019)
## 			Checked the accessions and 4 of 5 are from Malaysia, and 1 is from Thailand. 
## 				Given that Gekko lionotum is not listed as being in Malaysia or Thailand, and that G. cicakterbang is listed as being in S. Thailand and Malaysia but not in Myanmar. We probably have 4/5 lionotum loci and possibly 1/5 cicakterbang loci.
## 			SOLUTION: Should map to both
## 			NEW SOLUTION: Have it map only to Gekko lionotum
## 			Set Gekko cicakterbang to NA in tree

## Ptychozoon kuhli
## 		- Ptychozoon kuhli -> Gekko kuhli, Gekko nicobarensis
## 		Checked all accessions, and all were collected in Malaysia. Therefore, this is Gekko kuhli (nicobarensis is endemic to the Nicobar islands)
## 			SOLUTION: Remove from Gekko nicobarensis, can rename in tree from Ptychozoon to Gekko kuhli.

## Uma notata and Uma rufopunctata both map as: 
##	 	- Uma notata -> Uma notata
##		- Uma cowlesi -> Uma cowlesi
##	 	- Uma rufopunctata -> Uma notata
##		- Uma rufopunctata -> Uma cowlesi
##			SOLUTION: Let notata be notata, cowlesi be cowlesi, and add Uma rufopunctata as an additional taxon.


#########
# Single species matches: Flagged because these seem like taxa that would have been changed previously.
lapply(singleMatches, function(x) which(sapply(resList, function(y) x %in% y) == TRUE))



# ------------------------------------------------
# Now, we implement these manual adjustments

resList['Adelphicos newmanorum'] <- NA

resList['Anatololacerta anatolica'] <- 'Anatololacerta anatolica'
resList['Anatololacerta budaki'] <- 'Anatololacerta budaki'
resList['Anatololacerta danfordi'] <- 'Anatololacerta danfordi'
resList['Anatololacerta pelasgiana'] <- 'Anatololacerta pelasgiana'
new <- rep(NA, ncol(masterTable))
names(new) <- colnames(masterTable)
new['tree'] <- 'Anatololacerta oertzeni'
addons <- c(addons, list(new))		

resList['Anolis rodriguezii'] <- 'Anolis rodriguezii'

resList['Anolis breslini'] <- NA #not in tree

resList['Aspidoscelis sexlineatus'] <- 'Aspidoscelis sexlineatus'

resList['Cryptoblepharus plagiocephalus'] <- 'Cryptoblepharus plagiocephalus'
resList[['Cryptoblepharus australis']] <- NA # not in tree
resList['Cryptoblepharus pannosus'] <- 'Cryptoblepharus pannosus'
new <- rep(NA, ncol(masterTable))
names(new) <- colnames(masterTable)
new['tree'] <- 'Cryptoblepharus carnabyi'
addons <- c(addons, list(new))		

resList['Eumeces schneiderii'] <- 'Eumeces schneiderii'

resList['Hemidactylus triedrus'] <- 'Hemidactylus triedrus'
resList['Hemidactylus whitakeri'] <- 'Hemidactylus whitakeri'

resList['Hypsiglena unaocularus'] <- 'Hypsiglena unaocularus'
resList['Hypsiglena ochrorhynchus'] <- 'Hypsiglena ochrorhynchus'

resList['Lycodon subannulatus'] <- 'Lycodon subannulatus'
resList['Lycodon davisonii']

resList['Macrovipera lebetinus']
resList['Macrovipera schweizeri'] <- 'Macrovipera schweizeri'

resList['Madascincus polleni'] <- 'Madascincus polleni'
resList['Madascincus miafina'] <- NA # not in tree
new <- rep(NA, ncol(masterTable))
names(new) <- colnames(masterTable)
new['tree'] <- 'Madascincus intermedius'
addons <- c(addons, list(new))

resList['Magliophis exiguus'] <- 'Magliophis exiguus'

resList['Mochlus sundevallii'] <- 'Mochlus sundevallii'

resList['Sinonatrix percarinatus'] <- 'Sinonatrix percarinatus'

resList['Xenochrophis flavipunctatus']
resList['Xenochrophis piscator']
resList['Fowlea flavipunctatus']
resList['Fowlea piscator']
resList[['Fowlea melanzostus']] <- NA

resList['Gloydius intermedius'] <- 'Gloydius intermedius'
resList['Gloydius shedaoensis'] <- 'Gloydius shedaoensis'
new <- rep(NA, ncol(masterTable))
names(new) <- colnames(masterTable)
new['tree'] <- 'Gloydius saxatilis'
addons <- c(addons, list(new))		

resList['Kladirostratus acutus'] <- 'Kladirostratus acutus'
resList['Psammophylax acutus'] # should not be present
resList[['Kladirostratus togoensis']] <- NA

resList['Gekko cicakterbang'] <- NA


resList['Uma notata'] <- 'Uma notata'
resList[['Uma cowlesi']] <- NA # not present in tree 
resList['Uma rufopunctata'] # not present because we are adding this taxon
new <- rep(NA, ncol(masterTable))
names(new) <- colnames(masterTable)
new['tree'] <- 'Uma rufopunctata'
addons <- c(addons, list(new))




which(lengths(resList) > 1)
addons

dropFromTree <- c('Anolis rodriguezi', 'Aspidoscelis sexlineata', 'Eumeces schneideri', 'Hemidactylus subtriedrus', 'Hypsiglena ochrorhyncha', 'Magliophis exiguum', 'Mochlus sundevalli', 'Sinonatrix percarinata', 'Psammophylax acutus')

# check that these taxa to drop, are not listed elsewhere
for (i in 1:length(resList)) {
	resList[[i]] <- setdiff(resList[[i]], dropFromTree)
	if (length(resList[[i]]) == 0) {
		resList[[i]] <- NA
	}
}


# What changes to make to tree taxon names?
# Xenochrophis flavipunctatus renamed to Fowlea flavipunctatus
# Xenochrophis piscator renamed to Fowlea piscator
# Ptychozoon kuhli renamed to Gekko kuhli.
# Chilomeniscus stramineus renamed to Sonora straminea.

# a few more taxa that are being dropped from the tree
dropFromTree <- c('Asymblepharus sikimmensis', 'Correlophus belepensis', 'Eremias isfahanica', 'Leptophis modestus', 'Lophognathus gilberti', 'Pseudoxyrhopus analabe', 'Sibon noalamina')

# check that these taxa to drop, are not listed elsewhere
for (i in 1:length(resList)) {
	resList[[i]] <- setdiff(resList[[i]], dropFromTree)
	if (length(resList[[i]]) == 0) {
		resList[[i]] <- NA
	}
}

# Now there should only be a single taxon name mapped to any of the accepted Reptile DB names
# Fill in the table

identical(masterTable$repDB, names(resList))

masterTable$tree <- unlist(resList)

# add names not listed as accepted taxa
masterTable <- rbind.data.frame(masterTable, do.call(rbind, addons))

table(is.na(masterTable$tree))

# how many cases are there where the repDB name and tree name don't match?
table(apply(masterTable, 1, function(x) x[1] == x[2] | is.na(x[2])))
masterTable[!apply(masterTable, 1, function(x) x[1] == x[2] | is.na(x[2])),]

head(masterTable)


# useful code snippets
# which((sapply(resList, function(x) 'Adelphicos quadrivirgatus' %in% x)) == TRUE)
# taxonTable[which(taxonTable$repdbTaxon == tax), ]

# #############################################################
# ------------------------------------------------------------
# Carry out same process with geographic range taxa

resList <- replicate(length(repDB_acceptedBinomial), NA, simplify = FALSE)
names(resList) <- repDB_acceptedBinomial
addons <- list()

# we will also keep track of the mapping of Roll et al. taxa to geogTaxa as listed in the masterFile
RollConversion <- data.frame(Roll = geogTaxa, masterGeog = NA)

singleMatches <- c()
for (i in 1:length(geogTaxa)) {
	
	tax <- geogTaxa[i]
	
	#message('\t', tax)
	
	# is this taxon present in the reptile BD?
	if (tax %in% repDB_acceptedBinomial) {
		
		# we can treat this as a direct match
		ind <- which(masterTable$repDB == tax)
		resList[[ind]] <- c(setdiff(resList[[ind]], NA), tax)
	
	} else {
		
		# taxon is not immediately found in reptile DB
		
		# try matching
		match <- matchToReptileDB(tax, considerSynonyms = TRUE, resolveToSubSpecies = FALSE, returnMultipleMatches = TRUE)
		
		# was a match found via exact matching?
		if ('authority' %in% colnames(match)) {
			
			if (nrow(match) == 1) {
				# single match
				stop('single match accepted')
			}
			
			ind <- which(masterTable$repDB %in% match$taxon)
			for (j in ind) {
				resList[[j]] <- c(setdiff(resList[[j]], NA), tax)
			}
		}
		
		# was a match found via fuzzy matching and/or synonymy?
		if (!'authority' %in% colnames(match)) {				
			if (!all(is.na(match['repdbTaxon']))) {
				# something was found
				if (is.null(dim(match))) {
					# a single match was found
					singleMatches <- c(singleMatches, tax)
					# resList[[which(masterTable$repDB %in% match['repdbTaxon'])]] <- tax
					ind <- which(masterTable$repDB %in% match['repdbTaxon'])
					resList[[ind]] <- c(setdiff(resList[[ind]], NA), tax)
					
				} else if (!is.null(dim(match))) {
					# multiple matches found
					# in this case, we will fill in the table for multiple Rep DB taxa

					ind <- which(masterTable$repDB %in% match$taxon)
					for (j in ind) {
						resList[[j]] <- c(setdiff(resList[[j]], NA), tax)
					}
				}
			} else {
				# No matches found. We will add these to the table.
				addons
				new <- rep(NA, ncol(masterTable))
				names(new) <- colnames(masterTable)
				new['geog'] <- tax
				addons <- c(addons, list(new))		
			}
		}		
	}
}

which(lengths(resList) > 1)

# make table for manual eval
geogManual <- data.frame(RepDB = names(resList), geogTaxa = sapply(resList, function(x) paste0(x, collapse = ', ')), manualFix = '', geogUnion = '')
rownames(geogManual) <- NULL

# limit it to those rows with multiple options, or where repDB name does not match filled in name
checkInd <- c()
for (i in 1:length(resList)) {
	if (length(resList[[i]]) == 1) {
		if (!anyNA(resList[[i]]) & names(resList)[i] != resList[[i]]) {
			checkInd <- c(checkInd, i)
		}
	} else {
		#checkInd <- c(checkInd, i)
	}
}
spotCheck <- geogManual[checkInd,]

checkInd <- c()
for (i in 1:length(resList)) {
	if (length(resList[[i]]) == 1) {
		if (!anyNA(resList[[i]]) & names(resList)[i] != resList[[i]]) {
			#checkInd <- c(checkInd, i)
		}
	} else {
		checkInd <- c(checkInd, i)
	}
}

manual <- geogManual[checkInd,]


#write.csv(manual, 'geogTaxa_manual.csv')
#write.csv(spotCheck, 'geogTaxa_spotCheck.csv')

# new <- c()
# for (i in 1:nrow(manual1)) {
	# if (!(manual1[i,1] %in% manual[,1] & manual1[i,2] %in% manual[,2])) {
		# new <- c(new, i)
	# }
# }
# write.csv(manual1[new,], 'geogTaxa_manual_part2.csv')

which(lengths(resList) > 1) 
table(sapply(resList, function(x) x[1] == x[2] | is.na(x[2])))

# For geographic data, if 2 taxa have been synonymized, we may want to union their geographic ranges. Let's keep track here:
## Amphisbaena talisiae: A. mensae and A. talisiae. 
## 		These have been synonymized under A. talisiae. Should union ranges and remove A. mensae from geog names.

# Read in manually corrected dataset
manual <- read.csv('geogTaxa_manual_completed.csv')
manual2 <- read.csv('geogTaxa_manual_part2_completed.csv')
head(manual)
head(manual2)

table(lengths(resList[manual[,1]]))
table(lengths(resList[manual2[,1]]))

# each element in resList is an accepted RepDB name.
# the ones that were flagged for the manual fixes table were ones where more than one name matched
# We will replace resList[[ind]] with what was decided manually
# if any taxa were flagged to be geographically unioned, we will flag those again
geogUnion <- replicate(length(repDB_acceptedBinomial), NA, simplify = FALSE)
names(geogUnion) <- repDB_acceptedBinomial

for (i in 1:nrow(manual)) {
	
	ind <- which(names(resList) == manual[i, 'RepDB'])
	resList[ind]
	resList[[ind]] <- manual[i, 'RollTaxonFix']

	if (manual[i, 'geogUnion'] != '') {
		geogUnion[[ind]] <- manual[i, 'geogUnion']
	}		
}

for (i in 1:nrow(manual2)) {
	
	ind <- which(names(resList) == manual2[i, 'RepDB'])
	resList[ind]
	resList[[ind]] <- manual2[i, 'RollTaxonFix']

	if (manual2[i, 'geogUnion'] != '') {
		geogUnion[[ind]] <- manual2[i, 'geogUnion']
	}		
}



# deal with the taxa that were not matched to anything
do.call(rbind, addons)

# Afroablepharus seydeli
resList['Panaspis seydeli'] <- 'Afroablepharus seydeli'

# Ameivula littoralis
resList['Glaucomastix littoralis'] <- 'Ameivula littoralis'

# Amphisbaena maximus
resList['Leposternon maximus'] <- 'Amphisbaena maximus'

# Anolis ordinatus
## synonym of Anolis sagrei
# should geog union with A. sagrei
resList['Anolis sagrei']
geogUnion['Anolis sagrei'] <- 'Anolis sagrei, Anolis ordinatus'

# Anolis rhombifer
# synonym of Anolis lemurinus
resList['Anolis lemurinus']
geogUnion['Anolis lemurinus'] <- 'Anolis lemurinus, Anolis rhombifer'

# Anolis williamsii
# synonym of Anolis sericeus
resList['Anolis sericeus']
geogUnion['Anolis sericeus'] <- 'Anolis sericeus, Anolis williamsii'

# Basiliscus basilicus
# typo!
resList['Basiliscus basiliscus'] <- 'Basiliscus basilicus'

# Coluber barbouri
resList['Masticophis barbouri'] <- 'Coluber barbouri'

# Dryocalamus philippinus
resList['Lycodon philippinus'] <- 'Dryocalamus philippinus'

# Dryocalamus tristrigatus
resList['Lycodon tristrigatus'] <- 'Dryocalamus tristrigatus'

# Leiocephalus anonymous, Leiocephalus apertosulcus
resList['Leiocephalus carinatus']
geogUnion['Leiocephalus carinatus'] <- 'Leiocephalus carinatus, Leiocephalus anonymous, Leiocephalus apertosulcus'

# Leiocephalus partidus
resList['Leiocephalus partitus'] <- 'Leiocephalus partidus'

# Mochlus mabuiiforme
resList['Mochlus mabuiiformis'] <- 'Mochlus mabuiiforme'

# Orthriophis hodgsoni
resList['Elaphe hodgsoni'] <- 'Orthriophis hodgsoni'

# Panaspis nimbaensis
resList['Panaspis tristaoi'] <- 'Panaspis nimbaensis'

# Riama kiziriani
resList['Andinosaura kiziriani'] <- 'Riama kiziriani'

# Riama aurea
resList['Andinosaura aurea'] <- 'Riama aurea'

# Sphenomorphus aruensis
resList['Sphenomorphus jobiensis']
geogUnion['Sphenomorphus jobiensis'] <- 'Sphenomorphus jobiensis, Sphenomorphus aruensis'

# Sphenomorphus orientale
resList['Sphenomorphus orientalis'] <- 'Sphenomorphus orientale'

# Sphenomorphus striatopunctatum
resList['Lankascincus taprobanensis']
geogUnion['Lankascincus taprobanensis'] <- 'Lankascincus taprobanensis, Sphenomorphus striatopunctatum'

# Sphenomorphus textum
resList['Tytthoscincus textus'] <- 'Sphenomorphus textum'


# other manual adjustments -- singleGeog-to-multiRepDB is just problematic
resList['Andinosaura crypta'] <- 'Riama crypta'
resList['Diploderma drukdaypo'] <- NA
resList['Diploderma laeviventre'] <- NA
resList['Diploderma vela'] <- NA
resList['Dravidogecko douglasadamsi'] <- NA
resList['Dravidogecko janakiae'] <- NA
resList['Dravidogecko meghamalaiensis'] <- NA
resList['Dravidogecko septentrionalis'] <- NA
resList['Dravidogecko tholpalli'] <- NA
resList['Gekko cicakterbang'] <- NA
resList['Liopeltis pallidonuchalis'] <- NA
resList['Loxopholis osvaldoi'] <- NA
resList['Toxicocalamus ernstmayri'] <- NA
resList['Panaspis massaiensis'] <- NA
resList['Panaspis tsavoensis'] <- NA
resList['Pholidobolus paramuno'] <- NA
resList['Pholidobolus ulisesi'] <- NA
resList['Pholidoscelis umbratilis'] <- NA
resList['Smithophis atemporalis'] <- NA
resList['Sonora cincta'] <- NA
resList['Sonora fasciata'] <- NA
resList['Trachydactylus hajarensis'] <- NA





table(sapply(geogUnion, anyNA))
geogUnion <- geogUnion[which(sapply(geogUnion, anyNA) == FALSE)]

# double check: for geogUnion, the list names should be accepted names, but the names within the vectors should be from the Roll et al. dataset.
table(names(geogUnion) %in% names(resList))
table(setNames(unlist(strsplit(unlist(geogUnion), ', ')), NULL) %in% geogTaxa) # good

# Now there should only be a single taxon name mapped to any of the accepted Reptile DB names
# Fill in the table

identical(masterTable$repDB, names(resList)) # doesn't match because we added on a few from tree
table(names(resList) %in% masterTable$repDB) 

table(lengths(resList)) # we expect all of length 1 
table(sapply(resList, function(x) x[1] == x[2] | is.na(x[2])))

for (i in 1:length(resList)) {
	masterTable[which(masterTable$repDB == names(resList)[i]), 'geog'] <- resList[[i]]
}

masterTable$geog[which(masterTable$geog == '')] <- NA

# how many RepDB taxa don't have a geog range?
table(is.na(masterTable$geog))

# how many tree taxa have a range
sum(apply(masterTable[, c('tree', 'geog')], 1, function(x) !is.na(x[1]) & !is.na(x[2]))) / length(na.omit(masterTable$tree)) # 93%

# which taxa have tip in tree but no geog range?
nrow(masterTable[which(!is.na(masterTable$tree) & is.na(masterTable$geog)),])
head(masterTable[which(!is.na(masterTable$tree) & is.na(masterTable$geog)),], 10)
repDB_noGeog <- na.omit(masterTable[which(!is.na(masterTable$tree) & is.na(masterTable$geog)), 'repDB'])

# We will add the geogUnion taxa in as geogTaxa because that makes most sense. The taxa that correspond
# to the RepDB taxon.
names(geogUnion)
intersect(names(geogUnion), masterTable[grep(',', masterTable$geog), 1])
setdiff(names(geogUnion), masterTable[grep(',', masterTable$geog), 1])
for (i in 1:length(geogUnion)) {
	masterTable[which(masterTable[, 'repDB'] == names(geogUnion)[i]), 'geog'] <- geogUnion[[i]]
}


# we will deal with filling these geography blanks in a subsequent script. 

# autocountries <- vector('list', length(repDB_noGeog))
# for (i in 1:length(repDB_noGeog)) {
	
	# autocountries[[i]] <- getCountryFromSpecies(sp = repDB_noGeog[i], db = 'squam')

# }

# read in Roll et al. supplementary table that includes explanations for excluded taxa
rollSuppl <- as.data.frame(readxl::read_excel('41559_2017_332_MOESM3_ESM.xlsx'))
rollSuppl <- rollSuppl[which(!rollSuppl$Group %in% c('croc', 'turtle')),]
head(rollSuppl)

intersect(repDB_noGeog, rollSuppl$Binomial)
setdiff(rollSuppl$Binomial, repDB_noGeog)
setdiff(repDB_noGeog, rollSuppl$Binomial)

rollSuppl[which(rollSuppl$Binomial %in% intersect(repDB_noGeog, rollSuppl$Binomial)), ]

sort(rollSuppl[grep('Marine', rollSuppl[, 'Reason not mapped'], ignore.case = TRUE), 'Binomial'])

# which taxa in mastertable are missing range because those are explicitly not mapped?
intersect(masterTable[which(!is.na(masterTable$tree) & is.na(masterTable$geog)), 'repDB'], rollSuppl$Binomial)
# 68 taxa that are in tree but without geog do not have it because there were valid reasons. 

# write these to file
write.csv(masterTable[which(!is.na(masterTable$tree) & is.na(masterTable$geog)),], 'inTreeNotGeog_verify.csv')


head(masterTable)

# Geog dataset does not contain marine taxa. Will need to code those separately. 
# Should be able to fill in the missing tree taxa manually

# should also be able to fill in from the eco table Geographic range, which seems to already have pulled in Roll et al. data.

# is it possible to reconnect geogTaxa to Roll et al. taxon names?
table(na.omit(masterTable$geog) %in% geogTaxa)
setdiff(masterTable$geog, geogTaxa)
table(na.omit(unlist(strsplit(masterTable$geog, ', '))) %in% geogTaxa) # all taxa in geog column are found in Roll et al. names


# #############################################################
# ------------------------------------------------------------
# Carry out same process with ecological data taxa

resList <- replicate(length(repDB_acceptedBinomial), NA, simplify = FALSE)
names(resList) <- repDB_acceptedBinomial
addons <- list()

singleMatches <- c()
for (i in 1:length(ecoTaxa)) {
	
	tax <- ecoTaxa[i]
	
	#message('\t', tax)
	
	# is this taxon present in the reptile BD?
	if (tax %in% repDB_acceptedBinomial) {
		
		# we can treat this as a direct match
		ind <- which(masterTable$repDB == tax)
		resList[[ind]] <- c(setdiff(resList[[ind]], NA), tax)
	
	} else {
		
		# taxon is not immediately found in reptile DB
		
		# try matching
		match <- matchToReptileDB(tax, considerSynonyms = TRUE, resolveToSubSpecies = FALSE, returnMultipleMatches = TRUE)
		
		# was a match found via exact matching?
		if ('authority' %in% colnames(match)) {
			
			if (nrow(match) == 1) {
				# single match
				stop('single match accepted')
			}
			
			ind <- which(masterTable$repDB %in% match$taxon)
			for (j in ind) {
				resList[[j]] <- c(setdiff(resList[[j]], NA), tax)
			}
		}
		
		# was a match found via fuzzy matching and/or synonymy?
		if (!'authority' %in% colnames(match)) {				
			if (!all(is.na(match['repdbTaxon']))) {
				# something was found
				if (is.null(dim(match))) {
					# a single match was found
					singleMatches <- c(singleMatches, tax)
					ind <- which(masterTable$repDB %in% match['repdbTaxon'])
					# resList[[which(masterTable$repDB %in% match['repdbTaxon'])]] <- tax
					resList[[ind]] <- c(setdiff(resList[[ind]], NA), tax)
					
				} else if (!is.null(dim(match))) {
					# multiple matches found
					
					# in this case, we will fill in the table for multiple Rep DB taxa

					ind <- which(masterTable$repDB %in% match$taxon)
					for (j in ind) {
						resList[[j]] <- c(setdiff(resList[[j]], NA), tax)
					}
				}
			} else {
				# No matches found. We will add these to the table.
				addons
				new <- rep(NA, ncol(masterTable))
				names(new) <- colnames(masterTable)
				new['eco'] <- tax
				addons <- c(addons, list(new))		
			}
		}		
	}
}

# which could not be linked to a single repDB taxon?
which(lengths(resList) > 1)

# which couldn't be associated with anything?
do.call(rbind, addons)

# --------------------------
# manually deal with these

which((sapply(resList, function(x) 'Aspidoscelis sonorae' %in% x)) == TRUE)
resList['Amphisbaena talisiae']

# Amphisbaena talisiae
# has been synonymized with A. mensae. Will keep A. talisiae, and will drop mensae.
resList['Amphisbaena talisiae'] <- 'Amphisbaena talisiae'

# Anolis breslini, Anolis whitemani, Anolis saxatilis
## A. breslini and A. saxatilis are accepted RepDB taxa
## A. breslini was a subspecies of A. whitemani, but was elevated to species status 15 years ago
## A. whitemani has been synonymized with A. saxatilis and is no longer an accepted name.
## => Remove A. whitemani from A. breslini
resList['Anolis breslini'] <- 'Anolis breslini'


# Aspidoscelis sonorae
## A. flagellicaudus synonymized with A. sonorae
## As A. flagellicaudus not an accepted name and A. sonorae already in eco data, we drop A. flagellicaudus
resList['Aspidoscelis sonorae'] <- 'Aspidoscelis sonorae'

# Brookesia antakarana
## B. ambreensis has been synonymized with B. antakarana
resList['Brookesia antakarana'] <- 'Brookesia antakarana'

# Celestus warreni
## C. carraui is now a subspecies of Celestus warreni, so we drop C. carraui
resList['Celestus warreni'] <- 'Celestus warreni'

# Chondrodactylus fitzsimonsi
## Chondrodactylus fitzsimonsi was a subspecies of Pachydactylus laevigatus
## Pachydactylus laevigatus seems to now be Chondrodactylus laevigatus (already is in resList)
resList['Chondrodactylus fitzsimonsi'] <- 'Chondrodactylus fitzsimonsi'
resList['Chondrodactylus laevigatus']

# Chondrodactylus pulitzerae
resList['Chondrodactylus pulitzerae'] <- 'Chondrodactylus pulitzerae'

# Chondrodactylus turneri
resList['Chondrodactylus turneri'] <- 'Chondrodactylus turneri'

# Cnemaspis upendrai
## Cnemaspis clivicola and Cnemaspis upendrai have been synonymized. Drop C. clivicola.
resList['Cnemaspis upendrai'] <- 'Cnemaspis upendrai'

# Diploderma iadinum
resList['Diploderma iadinum'] <- 'Japalura iadina'

# Diploderma swild
## D. swild is a new species. We will just map Japalura dymondi to D. dymondi and Japalura flaviceps to D. flaviceps.
resList['Diploderma swild'] <- NA

# Diploderma swinhonis
resList['Diploderma swinhonis'] <- 'Japalura swinhonis'

# Eumeces blythianus
## Eumeces schneideri blythianus has been elevated to species status
resList['Eumeces blythianus'] <- 'Eumeces blythianus'

# Eutropis macrophthalma
## E. grandis has been synonymized with E. macropthalma. Drop E. grandis.
resList['Eutropis macrophthalma'] <- 'Eutropis macrophthalma'

# Gonatodes humeralis
## Cnemaspis timoriensis synonym of Gonatodes humeralis, but uncertain. Drop Cnemaspis timoriensis.
resList['Gonatodes humeralis'] <- 'Gonatodes humeralis'

# Lankascincus fallax
## L. deraniyagalae has been synonymized with L. fallax
resList['Lankascincus fallax'] <- 'Lankascincus fallax'

# Lankascincus sripadensis
resList['Lankascincus sripadensis'] <- 'Lankascincus sripadensis'

# Lankascincus taprobanensis
## A bit unclear. L. munindradasai seems to be synonymized with L. taprobanensis
resList['Lankascincus taprobanensis'] <- 'Lankascincus taprobanensis'

# Lygodactylus picturatus
## Lygodactylus luteopicturatus has been synonymized with Lygodactylus picturatus
resList['Lygodactylus picturatus'] <- 'Lygodactylus picturatus'

# Mochlus guineensis
resList['Mochlus guineensis'] <- 'Mochlus guineensis'

# Mochlus sundevallii
## Mochlus afer is synonymized with M. sundevallii
resList['Mochlus sundevallii'] <- 'Mochlus sundevalli'

# Phymaturus spurcus
## P. spectabilis has been synonymized with P. spurcus. Remove P. spectabilis
resList['Phymaturus spurcus'] <- 'Phymaturus spurcus'

# Sitana ponticeriana
## Sitana bahiri has been synonymized with S. ponticeriana
resList['Sitana ponticeriana'] <- 'Sitana ponticeriana'

# Acontias plumbeus
## Acontias poecilus synonymized with Acontias plumbeus. Remove A. poecilus. 
resList['Acontias plumbeus'] <- 'Acontias plumbeus'

# Agama agama
## Agama sylvana has been synonymized with A. africana. Since both A. agama and A. africana already present in eco dataset, we will drop A. sylvana. 
resList['Agama agama'] <- 'Agama agama'

# Amphisbaena darwinii
## A. heterozonata and A. trachura have been synonymized with A. darwinii.
resList['Amphisbaena darwinii'] <- 'Amphisbaena darwinii'

# Anolis cybotes
## A. haetianus has been synonymized with A. cybotes.
resList['Anolis cybotes'] <- 'Anolis cybotes'

# Anolis fuscoauratus
## A. scapularis synonymized with A. fuscoauratus
resList['Anolis fuscoauratus'] <- 'Anolis fuscoauratus'

# Cnemaspis goaensis
## C. indraneildasii synonymized with A. goaensis
resList['Cnemaspis goaensis'] <- 'Cnemaspis goaensis'

# Delma butleri
## Delma haroldi synonym of butleri
resList['Delma butleri'] <- 'Delma butleri'

# Diploderma laeviventre
## The taxonomy of Diploderma and Japalura is complicated. 
## Simplest would be to map J. flaviceps to D. flaviceps and J. laeviventris to D. laeviventre.
resList['Diploderma laeviventre'] <- 'Japalura laeviventris'
resList['Diploderma flaviceps'] <- 'Japalura flaviceps'

# Diploderma vela
## Should link Japalura vela to Diploderma vela
resList['Diploderma vela'] <- 'Japalura vela'

# Gekko nicobarensis
## Should link Ptychozoon nicobarensis to Gekko nicobarensis, and P. kuhli to Gekko kuhli
resList['Gekko nicobarensis'] <- 'Ptychozoon nicobarensis'

# Liolaemus audituvelatus
## L. manueli synonymized with audituvelatus
resList['Liolaemus audituvelatus'] <- 'Liolaemus audituvelatus'

# Liolaemus leopardinus
## L. ramonensis synonymized with leopardinus
resList['Liolaemus leopardinus'] <- 'Liolaemus leopardinus'

# Mochlus paedocarinatus
## Mochlus laeviceps paedocarinatum => now Mochlus paedocarinatus
## but Mochlus laeviceps still valid species
resList['Mochlus paedocarinatus'] <- 'Mochlus paedocarinatum'

# Pachydactylus atorquatus
## P. goodi synonymized with atorquatus
resList['Pachydactylus atorquatus'] <- 'Pachydactylus atorquatus'

# Tropidurus lagunablanca
## T. tarara and T. teyumirim synonymized with T. lagunablanca
resList['Tropidurus lagunablanca'] <- 'Tropidurus lagunablanca'

# Tytthoscincus butleri
## T. langkawiensis synonymized with T. butleri
resList['Tytthoscincus butleri'] <- 'Tytthoscincus butleri'


### The following were not matched to anything. We will verify.
do.call(rbind, addons)

# Leiocephalus anonymous, apertosulcus, partidus are extinct. 
# Leiocephalus partidus still in repDB so might as well correct.
# Ptyodactylus togensis is synonym of Ptyodactylus ragazzii
# Pholidoscelis turukaeraensis also extinct. 
resList['Leiocephalus partitus'] <- 'Leiocephalus partidus'
resList['Mochlus mabuiiformis'] <- 'Mochlus mabuiiforme'
resList['Panaspis tristaoi'] <- 'Panaspis nimbaensis'
resList['Phelsuma v-nigra'] <- 'Phelsuma vnigra'
resList['Sphenomorphus orientalis'] <- 'Sphenomorphus orientale'




# Now there should only be a single taxon name mapped to any of the accepted Reptile DB names
table(lengths(resList))
table(sapply(resList, function(x) x[1] == x[2] | is.na(x[2])))

# Fill in the table

identical(masterTable$repDB, names(resList)) # doesn't match because we added on a few from tree
table(names(resList) %in% masterTable$repDB) 

for (i in 1:length(resList)) {
	masterTable[which(masterTable$repDB == names(resList)[i]), 'eco'] <- resList[[i]]
}

# is it possible to reconnect ecoTaxa to Meiri 2018 taxon names?
table(na.omit(masterTable$eco) %in% ecoTaxa)
setdiff(na.omit(masterTable$eco), ecoTaxa)
table(na.omit(unlist(strsplit(masterTable$eco, ', '))) %in% ecoTaxa) # all taxa in eco column are found in Meiri 2018 names


###############################################
# of the species without geographic range data, how many have geographic range data in the eco dataset, in terms of country?
ind <- which(is.na(masterTable$geog) & !is.na(masterTable$eco))
eco[which(eco$Binomial %in% masterTable[ind, 'eco']), 'Geographic.Range']
# quite a few
###############################################
## Fill in family / genus

genusNotFound <- c()
for (i in 1:nrow(masterTable)) {
	
	tax <- masterTable[i, 'repDB']
	genus <- strsplit(tax, '\\s+|_')[[1]][1]
	masterTable[i, 'genus'] <- genus
	if (genus %in% taxonTable$genus) {
		fam <- taxonTable[which(taxonTable$genus == genus)[1], 'family']
		masterTable[i, 'family'] <- fam
	} else {
		genusNotFound <- c(genusNotFound, genus)	
	}
}

for (i in 1:nrow(masterTable)) {
	
	if (is.na(masterTable[i, 'repDB'])) {
		tax <- masterTable[i, 'tree']
		genus <- strsplit(tax, '\\s+|_')[[1]][1]
		masterTable[i, 'genus'] <- genus
		if (genus %in% taxonTable$genus) {
			fam <- taxonTable[which(taxonTable$genus == genus)[1], 'family']
			masterTable[i, 'family'] <- fam
		} else {
			genusNotFound <- c(genusNotFound, genus)	
		}
	}
}

table(is.na(masterTable$genus))
table(is.na(masterTable$family))

genusNotFound <- sort(unique(genusNotFound))

table(genusNotFound %in% pyronTaxonomy$Genus)

for (i in 1:length(genusNotFound)) {	
	if (genusNotFound[i] %in% pyronTaxonomy$Genus) {	
		fam <- pyronTaxonomy[which(pyronTaxonomy$Genus == genusNotFound[i])[1], 'Family']
		masterTable[which(masterTable$genus == genusNotFound[i]), 'family'] <- fam
	}
}

genusNotFound <- sort(unique(masterTable[which(is.na(masterTable$family)), 'genus']))

masterTable[which(masterTable$genus == 'Arcanumophis'), 'family'] <- 'Colubridae'
masterTable[which(masterTable$genus == 'Austroablepharus'), 'family'] <- 'Scincidae'
masterTable[which(masterTable$genus == 'Cenaspis'), 'family'] <- 'Colubridae'
masterTable[which(masterTable$genus == 'Dendrosauridion'), 'family'] <- 'Gymnophthalmidae'
masterTable[which(masterTable$genus == 'Dravidogecko'), 'family'] <- 'Gekkonidae' 
masterTable[which(masterTable$genus == 'Kuniesaurus'), 'family'] <- 'Scincidae'
masterTable[which(masterTable$genus == 'Microauris'), 'family'] <- 'Agamidae'
masterTable[which(masterTable$genus == 'Monilesaurus'), 'family'] <- 'Agamidae'
masterTable[which(masterTable$genus == 'Proahaetulla'), 'family'] <- 'Colubridae'
masterTable[which(masterTable$genus == 'Smithophis'), 'family'] <- 'Colubridae'

table(is.na(masterTable$genus))
table(is.na(masterTable$family))

for (i in 1:ncol(masterTable)) {
	masterTable[, i] <- gsub(',_', ', ', gsub('\\s+', '_', masterTable[, i]))
}

head(masterTable)
tail(masterTable)

write.csv(masterTable, outfile, row.names = FALSE)


# save geogUnion
saveRDS(geogUnion, '~/Dropbox/squamatePhylo/2019/empirical/geogUnion.rds')



# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------

# Read in metadata associated with phylogenetic inference

basedir <- '~/Dropbox/Oz_Crown_Ages/'
metaTableFile <- paste0(basedir, 'phylogenetic_inference/fasta_March2020/metatable3.csv')
taxonTable <- read.csv(metaTableFile, stringsAsFactors=FALSE)
head(taxonTable)

colnames(masterTable)[which(colnames(masterTable) == 'tree')] <- 'treename'

# add to our master taxonomy table whether taxa have genetic data, genomic data, or are imputed, and whether taxon is in topological constraint.
# A taxon can be in the topological constraint but be listed as "genbank" because it may be a taxon that was UCE or AHE only. 

treeDataType <- character(nrow(masterTable))
inTopoConstraint <- logical(nrow(masterTable))
for (i in 1:nrow(masterTable)) {
	
	sp <- masterTable[i, 'treename']
	
	if (!is.na(sp)) {
#		if (taxonTable[which(taxonTable$finalTipName == sp)[1], 'type'] == 'GenBank') {
		if (all(taxonTable[which(taxonTable$finalTipName == sp), 'type'] == 'GenBank')) {
			treeDataType[i] <- 'genbank'
		} else {
			treeDataType[i] <- 'genomic'
		}
		inTopoConstraint[i] <- taxonTable[which(taxonTable$finalTipName == sp)[1], 'inTopoConstraint']
	} else {
		treeDataType[i] <- 'imputed'
	} 
}

# for those that are in topological constraint but that used genbank data for the tree inference, let's change the data type flag to genomic, since those taxa really do have genomic data. 

#ind <- which(treeDataType == 'genbank' & inTopoConstraint == TRUE)
#treeDataType[ind] <- 'genomic'

table(treeDataType)
table(inTopoConstraint)

masterTable$dataType <- treeDataType
masterTable$inTopoConstraint <- inTopoConstraint


# add flag for "is one-to-one with tree" for correspondance between tree and reptileDB taxonomy. This will make it easier to recognize such issues without analyzing in R. 

masterTable[which(duplicated(masterTable[, 'treename'], incomparables = NA)), c('repDB', 'treename')]
masterTable[which(duplicated(masterTable[, 'repDB'], incomparables = NA)), c('repDB', 'treename')]

ind <- which(duplicated(masterTable[, 'treename'], incomparables = NA))
ind <- which(masterTable$treename %in% masterTable[ind, 'treename'])
masterTable$oneToOne <- TRUE
masterTable$oneToOne[ind] <- FALSE

# colnames(masterTable)[grep('tree', colnames(masterTable))] <- 'treename'

masterTable <- masterTable[order(masterTable$family, masterTable$genus, masterTable$repDB), ]

## Write core organizational file: 
### master table with RepDB, tree, genus, family, data type and constraint. 
### This one is organized around ReptileDB taxonomy as main column.
head(masterTable)

write.csv(masterTable, paste0(basedir, 'taxon-attributes/squamatesCoreTaxonomy_repDB.csv'), row.names = FALSE)


## Write core organizational file organized around tree taxa
# This one is similar, but we will fill in tree taxon NA's with repDB names as those are already flagged as imputed.

ind <- which(is.na(masterTable$treename))
table(masterTable[ind, 'dataType'])
table(masterTable[setdiff(1:nrow(masterTable), ind), 'dataType'])

masterTable[ind, 'treename'] <- masterTable[ind, 'repDB']
masterTable <- masterTable[, c('treename', 'family', 'genus', 'dataType', 'inTopoConstraint', 'oneToOne')]

anyDuplicated(masterTable)
anyDuplicated(masterTable$treename)

masterTable <- masterTable[!duplicated(masterTable$treename),]

# a few taxa have different genera in the treename compared to the official RepDB name.
## as this file is organized around the tree, have the genus match the treename
ind <- which(!(gsub('(.+)_(.+)', '\\1', masterTable$treename) == masterTable$genus))
masterTable[ind,]
masterTable$genus <- gsub('(.+)_(.+)', '\\1', masterTable$treename)

head(masterTable)
nrow(masterTable)

write.csv(masterTable, paste0(basedir, 'taxon-attributes/squamatesCoreTaxonomy_treeTaxa.csv'), row.names = FALSE)



