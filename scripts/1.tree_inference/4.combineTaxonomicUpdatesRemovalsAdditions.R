# Given automated 1-to-1 matching, manual corrections and adjustments, merge all of this information into a simpler table containing:
#	- accession
# 	- locus name
# 	- ncbi taxon
# 	- ncbi taxon id
# 	- updated (probably ReptileDB) taxon name
# 	- whether it is an addition, an update, or a seq slated for removal.


# this is more of a master file from a previous step. It contains all accessions currently in use. 
# We will make a few adjustments to this file first, and then we will extract the relevant updates.

dat <- read.csv('~/Dropbox/squamatePhylo/2019/topologyScripts/allAccessionsMatchingToReptileDB.csv', stringsAsFactors = FALSE)

head(dat)

options(warn=1)

#################
# No matches

# First, we will load in manual adjustments for ncbi taxa that were not matching to anything.
noMatchTaxa <- read.csv('~/Dropbox/Oz_Crown_Ages/ncbi_noMatchTaxa.csv', stringsAsFactors=FALSE)
noMatchTaxa <- unique(noMatchTaxa)
head(noMatchTaxa)

# incorporate this information into dat
for (i in 1:nrow(noMatchTaxa)) {
	
	cat(i, ' ')
	id <- noMatchTaxa[i, 'ncbiTaxonID']
	
	if (is.na(id)) {
		if (is.na(noMatchTaxa[i, 'ncbiTaxon_ssp']) | noMatchTaxa[i, 'ncbiTaxon_ssp'] == '') {
			ind <- which(dat$ncbiTaxon == noMatchTaxa[i, 'ncbiTaxon'])
		} else {
			ind <- which(dat$ncbiTaxon == noMatchTaxa[i, 'ncbiTaxon'] & dat$ncbiTaxon_ssp == noMatchTaxa[i, 'ncbiTaxon_ssp'])
		}
	} else {
		ind <- which(dat$ncbiTaxonID == id)
	}
	
	if (length(ind) == 0) stop()
	
	dat[ind, 'repdbTaxon'] <- noMatchTaxa[i, 'manualIdentification']
	dat[ind, 'reason'] <- 'manual adjustment'
}

table(dat$reason, useNA = 'always')


###################
# Tonini and Zaher

toniniZaher <- read.csv('~/Dropbox/Oz_Crown_Ages/toniniZaher_resolution.csv', stringsAsFactors=FALSE)
head(toniniZaher)
nrow(toniniZaher)

anyNA(toniniZaher$manualIdentification)

# accession + locus provide unique identifiers
nrow(dat) == nrow(unique(dat[, c('accession', 'locus')]))

for (i in 1:nrow(toniniZaher)) {
	
	cat(i, ' ')	
	ind <- which(dat$accession == toniniZaher[i, 'accession'] & dat$locus == toniniZaher[i, 'locus'])
	dat[ind, 'repdbTaxon'] <- toniniZaher[i, 'manualIdentification']

	comm <- toniniZaher[i, 'comments']
	if (comm == '') {
		comm <- 'manually evaluated due to overlap with Tonini or Zaher'
	} else {
		comm <- paste0('manually evaluated due to overlap with Tonini or Zaher; ', comm)
	}
	
	dat[ind, 'reason'] <- comm
	
}

table(dat$reason, useNA = 'always')



#############
# SqCL taxa
# not all taxa used in genomic analyses are internal SqCL, so we will only be able to update taxonomy for internal SqCL

sqcl <- read.csv('~/Dropbox/Oz_Crown_Ages/sqCL_taxa_repDBmatching.csv', stringsAsFactors = FALSE)
head(sqcl)
nrow(sqcl)
table(sqcl$remainsUncertain, useNA = 'always')
table(sqcl$manualMatch == '')

missingTaxa <- c()
for (i in 1:nrow(sqcl)) {
	
	if (sqcl[i, 'manualMatch'] != '') {
		
		cat(i, ' ')
		
		#ind <- intersect(grep('sqcl', dat$accession, ignore.case = TRUE), which(dat$ncbiTaxon == gsub('_', ' ', sqcl[i, 'species'])))
		ind <- intersect(grep('sqcl', dat$accession, ignore.case = TRUE), grep(sqcl[i, 'sample'], dat$accession, ignore.case = TRUE))
	
		if (length(ind) > 0) {
			
			# grep does not seem to be accidentally matching multiple samples
			if (nrow(unique(dat[ind, c('ncbiTaxon', 'ncbiTaxon_ssp')])) > 1 & !i %in% c(985, 1089, 1109, 1114)) stop()
			
			dat[ind, 'repdbTaxon'] <- sqcl[i, 'manualMatch']
			dat[ind, 'reason'] <- paste0('manually evaluated; ', sqcl[i, 'manualNotes'])
		}
		
		if (length(ind) == 0) {
			missingTaxa <- c(missingTaxa, i)
		}
	}
}

# Missing taxa are those that are in our phylogenomic dataset, but are not SqCL technically. They are already on genbank.

# sqcl[missingTaxa,]

##########################
# Lack of 1-to-1 matching

# For records that were not able to be 1-to-1 mapped to ReptileDB (due either to complete lack of matching, or due to multiple hits), we will keep the original NCBI taxon name. 
table(dat[which(is.na(dat$repdbTaxon) | dat$repdbTaxon == ''), 'reason'], useNA = 'always')


# Which ones have been examined but are still blank?
unique(dat[intersect(which(dat$reason == 'manual adjustment'), which(is.na(dat$repdbTaxon) | dat$repdbTaxon == '')), c('ncbiTaxon','ncbiTaxon_ssp')])

# these are all taxa that look like "Gelanesaurus sp", so acceptable.

ind <- which(is.na(dat$repdbTaxon) | dat$repdbTaxon == '')

dat[ind, 'repdbTaxon'] <- dat[ind, 'ncbiTaxon']

table(dat[which(is.na(dat$repdbTaxon) | dat$repdbTaxon == ''), 'reason'], useNA = 'always')

# A few additional manual edits
# change all Chionactis occipitalis to Sonora occipitalis
# change all Parias flavomaculatus to Trimeresurus flavomaculatus
# change all Corallus hortulanus to Corallus hortulana

manualTaxonUpdates <- rbind(
						c('Chionactis occipitalis', 'Sonora occipitalis'),
						c('Parias flavomaculatus', 'Trimeresurus flavomaculatus'),
						c('Corallus hortulanus', 'Corallus hortulana')
)

for (i in 1:nrow(manualTaxonUpdates)) {
	dat[which(dat$repdbTaxon == manualTaxonUpdates[i,1]), 'repdbTaxon'] <- manualTaxonUpdates[i,2]
}



###################################
# all above this point involved updates to existing fields in the table. Below this point involve additions or removals from the table.
###################################

###############################################
# OTHER ADDITIONS/CORRECTIONS

# Rabosky, D.L., Hutchinson, M.N., Donnellan, S.C., Talaba, A.L. & Lovette, I.J. (2014) Phylogenetic disassembly of species boundaries in a widespread group of Australian skinks (Scincidae: Ctenotus). Molecular Phylogenetics and Evolution, 77, 71–82.

# Many of the Ctenotus taxa revised in this paper are already updated in ReptileDB. For eastern and south-eastern localities, Ctenotus robustus should be reassigned to Ctenotus spaldingi. But this involves a level of complication such that we will not address this.

# CYTB accession for Ctenotus robustus is in agreement with Rabosky 2014. Same is true for accession for Ctenotus spaldingi. We are not using locus ATPSB. Therefore no modifications needed.

# Ctenotus superciliaris: adding a accession below


# Liu, S. & Rao, D. (2019) A new species of the genus Acanthosaura from Yunnan, China (Squamata, Agamidae). ZooKeys, 888, 105–132.
# https://zookeys.pensoft.net/article/38491/
# Acanthosaura crucigera CYTB should be AY572888, and AY572887 should be dropped, because it is no longer associated with that taxon. We will add/remove these below.

### Are there any accessions that are duplicated with the same gene?
table(duplicated(dat[,c('locus','accession')])) 


###################################
# adjustments from Zaher blacklist

library(AnnotationBustR)
library(reutils)
library(Biostrings)


# This table lists accessions shared by our dataset and Zaher's blacklist. Listed accessions should be dropped, and replaced with replacement accession when available. 
zaherBlacklist <- read.csv('~/Dropbox/squamatePhylo/2019/topologyScripts/ZaherBlacklisted_updates.csv', stringsAsFactors=FALSE)
head(zaherBlacklist)

newSeqs <- matrix(nrow = 0, ncol = ncol(dat)+1)
colnames(newSeqs) <- c(colnames(dat), 'seq')

# update info in table, identify seqs to remove, add new seqs to newSeqs list



accessionsToAdd <- zaherBlacklist$replacementAccession
accessionsToAdd <- accessionsToAdd[accessionsToAdd != '']

## add a few more accessions here

# Ctenotus superciliaris CYTB from Rabosky 2014
accessionsToAdd <- c(accessionsToAdd, 'KJ505688') 

# Acanthosaura crucigera CYTB from Liu and Rao 2019
accessionsToAdd <- c(accessionsToAdd, 'AY572888') 

sqclFiles <- list.files('~/Dropbox/squamatePhylo/2019/SqCL/toAdd_toDB', pattern = '\\.fa(sta)?', full.names = TRUE)

for (i in 1:length(accessionsToAdd)) {
	
	message('\t', i)
	
	newRec <- rep(NA, ncol(newSeqs))
	names(newRec) <- colnames(newSeqs)
	
	if (grepl('sqcl', accessionsToAdd[i], ignore.case = TRUE)) {
		
		# identify appropriate gene region
		ind <- which(sapply(gsub('\\.fa$', '', basename(sqclFiles)), function(x) grepl(x, accessionsToAdd[i])) == TRUE)
		gene <- names(ind)
		
		tmp <- readDNAStringSet(sqclFiles[ind])
		
		# pull out the specific seq
		newRec['filename'] <- sqclFiles[ind]
		ind <- grep(gsub(paste0('_', gene), '', accessionsToAdd[i]), names(tmp))
		qq <- names(tmp)[ind]
		taxon <- strsplit(strsplit(qq, '@')[[1]][2], ' ')[[1]][-1]
		newRec['ncbiTaxonID'] <- strsplit(qq, '@')[[1]][1]
		newRec[c('ncbiTaxon', 'repdbTaxon')] <- paste(taxon[1:2], collapse = ' ')
		newRec[c('ncbiTaxon_ssp', 'repdbTaxon_ssp')] <- taxon[3]
		newRec['locus'] <- gene
		newRec['accession'] <- strsplit(strsplit(qq, '@')[[1]][2], ' ')[[1]][1]
		newRec['reason'] <- 'exact mapping via accepted'
		newRec['seq'] <- as.character(tmp[[ind]])
		
	} else {
		
		Sys.sleep(0.5)
		acc <- efetch(uid = accessionsToAdd[i], db = 'nuccore', rettype = "fasta", retmode = "text")
		
		acc <- strsplit(content(acc), split = '\n\n')[[1]]
		acc <- strsplit(acc, '\n')[[1]]
		accName <- acc[1]
		acc <- acc[-1]
		acc <- paste(acc, collapse = '')
		
		ind <- which(zaherBlacklist$replacementAccession == accessionsToAdd[i])
		
		if (length(ind) > 0) {
			
			newRec[c('ncbiTaxon', 'repdbTaxon')] <- zaherBlacklist[ind, 'taxon']
			newRec['ncbiTaxonID'] <- zaherBlacklist[ind, 'taxonID']
			newRec['locus'] <- gsub('_cluster[0-9]+$', '', zaherBlacklist[ind, 'cluster'])
			newRec['accession'] <- accessionsToAdd[i]
			newRec['reason'] <- 'exact mapping via accepted'
			newRec['seq'] <- acc
			newRec['filename'] <- 'downloaded'
			
		} else if (accessionsToAdd[i] == 'KJ505688') {
			
			newRec[c('ncbiTaxon', 'repdbTaxon')] <- 'Ctenotus superciliaris'
			newRec['locus'] <- 'CYTB'
			newRec['accession'] <- 'KJ505688'
			newRec['reason'] <- 'exact mapping via accepted'
			newRec['seq'] <- acc
			newRec['filename'] <- 'downloaded'

		} else if (accessionsToAdd[i] == 'AY572888') {
			
			newRec[c('ncbiTaxon', 'repdbTaxon')] <- 'Acanthosaura crucigera'
			newRec['ncbiTaxonID'] <- '103692'
			newRec['locus'] <- 'CYTB'
			newRec['accession'] <- 'AY572888'
			newRec['reason'] <- 'exact mapping via accepted'
			newRec['seq'] <- acc
			newRec['filename'] <- 'downloaded'
		} else {
			stop()
		}
		
	}
	
	
	newSeqs <- rbind(newSeqs, newRec)
	
}



dat <- rbind.data.frame(dat, newSeqs[, 1:ncol(dat)])



## DROP BLACKLISTED ACCESSIONS FROM MASTER TABLE

toDrop <- zaherBlacklist$accession
toDrop <- toDrop[toDrop != '']

# Also drop Acanthosaura crucigera CYTB accession
toDrop <- c(toDrop, 'AY572887')

toDrop <- sapply(toDrop, function(x) which(dat$accession == x), USE.NAMES=FALSE)

dat <- dat[-toDrop, ]


## DROP TAXA THAT WILL NOT BE USEFUL (i.e., lacking species)

table(grep('\\ssp\\.?$', dat$repdbTaxon, value=TRUE))
setdiff(1:nrow(dat), grep('\\s', dat$repdbTaxon))

dat <- dat[- grep('\\ssp\\.?$', dat$repdbTaxon, value=FALSE),]

# There is a non-usable taxon named "Iranian colubrine"
table(grep('colubrine', dat$repdbTaxon, value=TRUE))
dat <- dat[- grep('colubrine', dat$repdbTaxon, value=FALSE),]

table(grep('\\ssp\\.?$', dat$repdbTaxon, value=TRUE))



## Apparently, there are still some seqs < 300 bp. We will pull up each sequence to get its length.

dat$seqLength <- NA

uniqueFiles <- unique(dat$filename)
uniqueFiles <- uniqueFiles[!grepl('downloaded', uniqueFiles)]
uniqueFileList <- vector('list', length(uniqueFiles))
names(uniqueFileList) <- uniqueFiles
for (j in 1:length(uniqueFiles)) {
	
	f <- readDNAStringSet(uniqueFiles[j])
	# rename seqs
	seqNames1 <- setdiff(names(f), grep('@', names(f), value = TRUE))
	seqNames2 <- grep('@', names(f), value = TRUE)
	if (length(seqNames1) > 0) {
		seqNames1Acc <- strsplit(seqNames1, '\\s+')
		seqNames1Acc <- sapply(seqNames1Acc, function(x) x[1])
		seqNames1Acc <- gsub('\\.[0-9]$', '', seqNames1Acc)
		names(f)[sapply(seqNames1, function(x) which(names(f) == x))] <- seqNames1Acc
	}
	if (length(seqNames2) > 0) {			
		qq <- strsplit(seqNames2, '@') # first part is NCBI taxonID, second part is sample name/locus, and taxon name (maybe)
		seqNames2Acc <- unlist(lapply(qq, function(x) strsplit(x[2], '\\s+')[[1]][1]))
		names(f)[sapply(seqNames2, function(x) which(names(f) == x))] <- seqNames2Acc
	}
	
	uniqueFileList[[j]] <- f
}
	
# for each record, identify the appropriate sequence data file, and pull out the relevant sequence
for (j in 1:nrow(dat)) {
	cat(j, ' ')
	rec <- dat[j,]
	ind <- which(names(uniqueFileList) == rec$filename)
	if (length(ind) == 0 & !grepl('download', rec$filename)) stop()
	if (grepl('download', rec$filename)) {
		seq <- newSeqs[which(newSeqs[, 'accession'] == rec$accession & newSeqs[, 'locus'] == rec$locus), 'seq']
	} else {
		if (rec$accession %in% names(uniqueFileList[[ind]])) {
			seq <- uniqueFileList[[ind]][rec$accession,]
		} else if (paste0(rec$accession, '_', rec$locus) %in% names(uniqueFileList[[ind]])) {
			seq <- uniqueFileList[[ind]][paste0(rec$accession, '_', rec$locus),]
		} else {
			stop()
		}
	}
	seq <- as.character(seq)
	dat[j, 'seqLength'] <- nchar(seq)
}

anyNA(dat$seqLength)
range(dat$seqLength)
ind <- which(dat$seqLength < 300)

# drop those seqs that are under 300 base pairs in length
dat <- dat[- ind, ]

rm(uniqueFiles, uniqueFileList)

###############################################
## Updates, additions and corrections now made



##########################################################################
## Next step: Reorganization

# Organize into list of tables, where each list element is new taxon name
# We are ignoring subspecies

# first, add a column that we will use below for flagging records to be kept or removed
dat$keep <- TRUE

dat_repDB <- split(dat, dat$repdbTaxon)


# Here multiple seqs for the same gene region might be present for a taxon because 
# - I included all potential pyphlawd clusters, so there could be more than one for particular genes
# - I was treating subspecies separately before, whereas now we are lumping subspecies into species

# How to pick one representative per gene per taxon when multiple are present
## - (1) if one accession was used by Tonini or Zaher or SqCL, keep that one
## - (2) if one is associated with a taxon that was precisely matched, keep that one.
## - (3) check seq length. If one is truly longer than the other, keep the longer one.

for (i in 1:length(dat_repDB)) {
	if (any(duplicated(dat_repDB[[i]]$locus))) {
		
		message('\t', i)
		
	#	if (!i %in% c(4)) stop()
		
		uniqueLoci <- unique(dat_repDB[[i]]$locus)
		table(dat_repDB[[i]]$locus)
		for (j in 1:length(uniqueLoci)) {
			
			ind <- which(dat_repDB[[i]]$locus == uniqueLoci[j])
			if (length(ind) > 1) {
				
#				if (length(ind) > 2) stop()

				# this locus appears more than once
				dups <- dat_repDB[[i]][ind, ]

				# get sequence lengths
				seqL <- numeric(nrow(dups))
				for (k in 1:nrow(dups)) {
					f <- readDNAStringSet(dups[k, 'filename'])

					# rename seqs
					seqNames1 <- setdiff(names(f), grep('@', names(f), value = TRUE))
					seqNames2 <- grep('@', names(f), value = TRUE)
					seqNames1Acc <- strsplit(seqNames1, '\\s+')
					seqNames1Acc <- sapply(seqNames1Acc, function(x) x[1])
					seqNames1Acc <- gsub('\\.[0-9]$', '', seqNames1Acc)
					names(f)[sapply(seqNames1, function(x) which(names(f) == x))] <- seqNames1Acc

					if (length(seqNames2) > 0) {
						
						qq <- strsplit(seqNames2, '@') # first part is NCBI taxonID, second part is sample name/locus, and taxon name (maybe)
						seqNames2Acc <- unlist(lapply(qq, function(x) strsplit(x[2], '\\s+')[[1]][1]))
						names(f)[sapply(seqNames2, function(x) which(names(f) == x))] <- seqNames2Acc
					}
										
					seqAcc <- dups[k, 'accession']
					if (!seqAcc %in% names(f)) {
						if (paste0(seqAcc, '_', dups[k, 'locus']) %in% names(f)) {
							seqAcc <- paste0(seqAcc, '_', dups[k, 'locus'])
						} else if (paste0(seqAcc, '.1') %in% names(f)) {
							seqAcc <- paste0(seqAcc, '.1')
						} else {
							stop()
						}
					}
					seqL[k] <- nchar(f[seqAcc, ])
				}
				
										
				# is one SqCL or used by Tonini/Zaher?
				isSqCL <- grepl('sqcl', dups$accession, ignore.case = TRUE)
				
				# was one used by Tonini or Zaher?
				usedToniniZaher <- !is.na(dups$toniniTaxon) | !is.na(dups$zaherTaxon)
				
				# was one mapped more precisely to a reptileDB taxon, or was one manually adjusted?
				taxonMapping <- grepl('exact|manual', dups$reason)
				
				if (sum(isSqCL) > 0) {
					keep <- which(isSqCL == TRUE)
				} else if (sum(usedToniniZaher) > 0) {
					keep <- which(usedToniniZaher == TRUE)
				} else if (sum(taxonMapping) > 0) {
					keep <- which(taxonMapping == TRUE)
				} else {
					keep <- 1:nrow(dups)
				}
				
				# We don't necessarily want to favor SqCL for mitochondrial genes
				# So remove this priority for mt genes
				if (dups[1, 'locus'] %in% c('COI', 'CYTB', 'ND1', 'ND2', 'ND4', 'rRNA_12S', 'rRNA_16S')) {
					if (sum(usedToniniZaher) > 0) {
						keep <- which(usedToniniZaher == TRUE)
					} else if (sum(taxonMapping) > 0) {
						keep <- which(taxonMapping == TRUE)
					} else {
						keep <- 1:nrow(dups)
					}
				}
				
				
				if (length(keep) > 1) {
					
						
					if (length(which(seqL[keep] == max(seqL[keep]))) == 1) {
						keep <- keep[which.max(seqL[keep])]
					} else if (length(which(seqL[keep] == max(seqL[keep]))) < length(keep)) {
						keep <- keep[which(seqL[keep] == max(seqL[keep]))]
					}

				
					if (length(keep) > 1) {
						# At this stage, seqs have the same lengths. 
						# Let's consider multiple criteria. 
						
						keep2 <- which(isSqCL == TRUE | usedToniniZaher == TRUE | taxonMapping == TRUE)
						
						# if none of the criteria are satisfied, then just evaluate all for length
						if (length(keep2) == 0) {
							keep2 <- 1:nrow(dups)
						} 
						
						# try again filtering by length
						if (length(which(seqL[keep2] == max(seqL[keep2]))) == 1) {
							keep2 <- keep2[which.max(seqL[keep2])]
						} else if (length(which(seqL[keep2] == max(seqL[keep2]))) < length(keep2)) {
							keep2 <- keep2[which(seqL[keep2] == max(seqL[keep2]))]
						}
						
						if (length(keep2) > 1) {
							# At this point, we can just pick one randomly
							keep <- keep2[1]
						} else {
							keep <- keep2
						}
					}
				}
					
				# flag all but selected record for removal
				dat_repDB[[i]][setdiff(ind, ind[keep]), 'keep'] <- FALSE
				
			}
		}
	}
}

sum(sapply(dat_repDB, nrow))

# Now we drop those records that are no longer needed.
for (i in 1:length(dat_repDB)) {
	if (any(dat_repDB[[i]]$keep == FALSE)) {
		if (all(dat_repDB[[i]]$keep == FALSE)) stop()
		dat_repDB[[i]] <- dat_repDB[[i]][which(dat_repDB[[i]]$keep == TRUE),]	
	}
}

sum(sapply(dat_repDB, nrow))


###################################
# Can we fill in any gaps?

# Using the accessions table from pyphlawd, we can look to see if any accessions exist on genbank that could fill some gaps in our dataset

taxonTableFile <- '~/Dropbox/Oz_Crown_Ages/Lepidosauria_8504.table'
taxonTable <- data.table::fread(taxonTableFile, data.table=FALSE)
colnames(taxonTable) <- c('V1', 'taxonID', 'locus', 'accession', 'taxon', 'definition', 'title')
taxonTable$repdbTaxon <- NA
head(taxonTable)

# First we need to match ncbi taxon names to ReptileDB taxa
# We will also apply the same no-match corrections from above.
# Again, we will keep NCBI taxon names for non-matched taxa.

uniqueTableTaxa <- taxonTable$taxon
uniqueTableTaxa <- gsub('_', ' ', uniqueTableTaxa)
uniqueTableTaxa <- gsub('\\s+', ' ', uniqueTableTaxa)
uniqueTableTaxa <- gsub('(^\\s)|(\\s$)', '', uniqueTableTaxa)
taxonTable$taxon <- uniqueTableTaxa
uniqueTableTaxa <- sort(unique(uniqueTableTaxa))
length(uniqueTableTaxa)

repDBdir <- '~/Dropbox/Oz_Crown_Ages/reptileDB/'

# Read in Reptile Database data
acceptedNamesTable <- read.csv(paste0(repDBdir, 'acceptedTaxonTable.csv'), stringsAsFactors = FALSE)
subspeciesTable <- read.csv(paste0(repDBdir, 'subspeciesTable.csv'), stringsAsFactors = FALSE)
synonymTable <- read.csv(paste0(repDBdir, 'synonymTable.csv'), stringsAsFactors = FALSE)

# limit synonyms to post-1980
synonymTable <- synonymTable[which(synonymTable$year >= 1980 | is.na(synonymTable$year)),]

source('~/Dropbox/Oz_Crown_Ages/reptileDB/matchToReptileDB.R')


for (i in 1:length(uniqueTableTaxa)) {
	
	message('\t', i, '\t', uniqueTableTaxa[i])

	ind <- which(taxonTable$taxon == uniqueTableTaxa[i])

	match <- matchToReptileDB(uniqueTableTaxa[i], considerSynonyms = TRUE, resolveToSubSpecies = FALSE)
	
	# if no match, but in table of taxa that we handled mannually, then use that
	if (is.na(match['repdbTaxon'])) {
		if (taxonTable[ind[1], 'taxonID'] %in% noMatchTaxa$ncbiTaxonID) {
			match['repdbTaxon'] <- noMatchTaxa[which(noMatchTaxa$ncbiTaxonID == taxonTable[ind[1], 'taxonID']), 'manualIdentification']
		} 
	}
	
	taxonTable[ind, 'repdbTaxon'] <- match['repdbTaxon']
}


# for any taxa that did not match, copy over NCBI taxon name
table(is.na(taxonTable$repdbTaxon))
taxonTable[which(is.na(taxonTable$repdbTaxon)), 'repdbTaxon'] <- taxonTable[which(is.na(taxonTable$repdbTaxon)), 'taxon']

# apply those same taxonomic updates
for (i in 1:nrow(manualTaxonUpdates)) {
	taxonTable[which(taxonTable$repdbTaxon == manualTaxonUpdates[i,1]), 'repdbTaxon'] <- manualTaxonUpdates[i,2]
}


# To avoid confusion, we will drop the Zaher blacklisted accessions. 
blacklisted <- zaherBlacklist$accession
blacklisted <- blacklisted[blacklisted != '']

# Also drop Acanthosaura crucigera CYTB accession
blacklisted <- c(blacklisted, 'AY572887')

table(blacklisted %in% taxonTable$accession) # all are present

taxonTable <- taxonTable[- which(taxonTable$accession %in% blacklisted),]

table(blacklisted %in% taxonTable$accession) # all are present
# removed.



# For each taxon, we will try to find accessions in the taxonTable that are currently missing from our dataset

# read in Dan Portik's table of gene names and abbreviations
geneTable <- data.table::fread('~/SuperCRUNCH/data/locus-search-terms/Locus-Search-Terms-Squamate-Markers.txt', header=FALSE, data.table=FALSE)
geneList <- list()
for (i in 1:nrow(geneTable)) {
	geneList[[i]] <- c(strsplit(geneTable[i,2], ';')[[1]], strsplit(geneTable[i,3], ';')[[1]])
}
names(geneList) <- geneTable[,1]
geneList <- lapply(geneList, function(x) setdiff(x, 'mitochondrion, complete genome'))
names(geneList)[which(names(geneList) == 'CMOS')] <- 'cmos'
names(geneList)[which(names(geneList) == 'CO1')] <- 'COI'
names(geneList)[which(names(geneList) == 'ECEL1')] <- 'ECEL'
names(geneList)[which(names(geneList) == 'LZTS1')] <- 'LZTSS1'
names(geneList)[which(names(geneList) == 'NT3')] <- 'NTF-3'
names(geneList)[which(names(geneList) == '12S')] <- 'rRNA_12S'
names(geneList)[which(names(geneList) == '16S')] <- 'rRNA_16S'
names(geneList)[which(names(geneList) == 'SNCAIP')] <- 'SINCAIP'
names(geneList)[which(names(geneList) == 'VCPIP')] <- 'VCPIP1'
for (i in 1:length(geneList)) {
	geneList[[i]] <- unique(c(names(geneList)[i], geneList[[i]]))
}

# since I have added the AHE locus names to the descriptions, we might as well include those too
ind <- which(names(geneList) == 'VCPIP1')
geneList[[ind]] <- c('AHE-L203', geneList[[ind]])
ind <- which(names(geneList) %in% c('DLL', 'DLL1'))
geneList[[ind]] <- c('AHE-L183', geneList[[ind]])
ind <- which(names(geneList) == 'MSH6')
geneList[[ind]] <- c('AHE-L50', geneList[[ind]])
ind <- which(names(geneList) == 'ADNP')
geneList[[ind]] <- c('AHE-L381', geneList[[ind]])
ind <- which(names(geneList) == 'ZEB2')
geneList[[ind]] <- c('AHE-L60', geneList[[ind]])
ind <- which(names(geneList) == 'CAND1')
geneList[[ind]] <- c('AHE-L288', geneList[[ind]])
ind <- which(names(geneList) == 'GHSR')
geneList[[ind]] <- c('AHE-L113', geneList[[ind]])
ind <- which(names(geneList) == 'RAG1')
geneList[[ind]] <- c('AHE-L259', geneList[[ind]])
ind <- which(names(geneList) == 'SLC8A3')
geneList[[ind]] <- c('AHE-L309', geneList[[ind]])

geneList <- lapply(geneList, function(x) setdiff(x, c('mitochondrion complete genome', 'mitochondrial DNA complete genome', 'mitochondrion partial genome')))

# Let's avoid the ND1 vs CAND1 hell by using a negative lookaround (?<!CA)ND1 instead of just ND1 (need to remember to set perl=TRUE)
geneList[['ND1']] <- gsub('ND1', '(?<!CA)ND1', geneList[['ND1']])

# modify for regex
geneList <- lapply(geneList, function(x) gsub('\\s+', '\\\\\\s', x))

# add bounds
geneList <- lapply(geneList, function(x) paste0('\\b', x, '\\b'))

geneList <- lapply(geneList, function(x) paste0(x, collapse = '|'))

# limit gene table to those loci we are looking for
all(sort(unique(dat$locus)) %in% names(geneList))
geneList <- geneList[names(geneList) %in% unique(dat$locus)]


addedAccessions <- matrix(nrow = 0, ncol = 4)
colnames(addedAccessions) <- c('repdbTaxon', 'locus', 'accession', 'seqLength')


for (i in 1:length(dat_repDB)) {
	
	dat_repDB[[i]]$locus
	missingLoci <- setdiff(names(geneList), dat_repDB[[i]]$locus)
	
	taxonTableSub <- taxonTable[which(taxonTable$repdbTaxon == dat_repDB[[i]]$repdbTaxon[1]),]
	
	if (length(missingLoci) > 0) {
		for (j in 1:length(missingLoci)) {
			
			ind <- grep(geneList[[missingLoci[j]]], taxonTableSub$definition, perl = TRUE, ignore.case = TRUE)
			# identified accessions may already be in use
			ind <- ind[!taxonTableSub[ind, 'accession'] %in% dat_repDB[[i]]$accession]
			
			if (length(ind) > 0) {
				message('\t', i, '\tnew accessions found for ', missingLoci[j])
				
				# if multiple, which is longest sequence?
				if (length(ind) > 1) {
					Sys.sleep(0.5)
					newAcc <- try(FindLongestSeq(taxonTableSub[ind, 'accession']))
					if (inherits(newAcc, 'try-error')) {
						attempts <- 0
						while (inherits(newAcc, 'try-error') & attempts < 100) {
							message('\t\treattempting...')
							newAcc <- try(FindLongestSeq(taxonTableSub[ind, 'accession']))
						}
						if (inherits(newAcc, 'try-error')) stop()
					}
					if (nrow(newAcc) > 1) {
						newAcc <- newAcc[which.max(newAcc$Length),]
					}
					newAcc <- newAcc$Accession
				} else {
					newAcc <- taxonTableSub[ind, 'accession']
				}
				if (length(newAcc) > 1) stop()
				
				addedAccessions <- rbind(addedAccessions, c(repdbTaxon = dat_repDB[[i]]$repdbTaxon[1], locus = missingLoci[j], accession = newAcc, seqLength = NA))
			}
		}
	}
}

nrow(addedAccessions)

# Acquire sequence data and check length. 

addedAccessionsSeq <- vector('list', nrow(addedAccessions))
for (i in 1:nrow(addedAccessions)) {
	cat(i, ' ')
	Sys.sleep(0.5)
	zz <- try(efetch(uid = addedAccessions[i, 'accession'], db = 'nuccore', rettype = "fasta", retmode = "text"))
	if (inherits(zz, 'try-error')) {
		attempts <- 0
		while (inherits(zz, 'try-error') & attempts < 10) {
			message('\t\treattempting...')
			zz <- try(efetch(uid = addedAccessions[i, 'accession'], db = 'nuccore', rettype = "fasta", retmode = "text"))
		}
		if (inherits(zz, 'try-error')) stop()

	}
	zz <- strsplit(content(zz), split = '\n\n')[[1]]
	zz <- strsplit(zz, '\n')[[1]]
	zzName <- zz[1]
	zz <- zz[-1]
	zz <- paste(zz, collapse = '')
	addedAccessionsSeq[[i]] <- zz
	names(addedAccessionsSeq)[i] <- zzName
}


addedAccessions[,'seqLength'] <- sapply(addedAccessionsSeq, nchar)

addedAccessions <- as.data.frame(addedAccessions)
addedAccessions$repdbTaxon <- as.character(addedAccessions$repdbTaxon)
addedAccessions$locus <- as.character(addedAccessions$locus)
addedAccessions$accession <- as.character(addedAccessions$accession)
addedAccessions$seqLength <- as.numeric(as.character(addedAccessions$seqLength))

head(addedAccessions)

# Drop seqs that are under 300 bp
range(addedAccessions$seqLength)
ind <- which(addedAccessions$seqLength >= 300)
addedAccessions <- addedAccessions[ind, ]
addedAccessionsSeq <- addedAccessionsSeq[ind]

# Some of these accessions may in fact already be in the dataset but with an updated taxon name.
alreadyIncluded <- c()
for (i in 1:nrow(addedAccessions)) {
	cat(i, ' ')
	tmp <- sapply(dat_repDB, function(x) any(x$accession == addedAccessions[i, 'accession'] & x$locus == addedAccessions[i, 'locus']))
	if (any(tmp)) {
		alreadyIncluded <- c(alreadyIncluded, i)
	}
}

addedAccessions <- addedAccessions[- alreadyIncluded, ]
addedAccessionsSeq <- addedAccessionsSeq[- alreadyIncluded]



# add these records to the main data table
for (i in 1:nrow(addedAccessions)) {
	
	cat(i, ' ')
	spInd <- which(names(dat_repDB) == addedAccessions[i, 'repdbTaxon'])
	
	newRec <- rep(NA, ncol(dat_repDB[[spInd]]))
	names(newRec) <- colnames(dat_repDB[[spInd]])

	tax <- taxonTable[which(taxonTable$accession == addedAccessions[i, 'accession']), 'taxon']
	tax <- strsplit(tax, split = '\\s+')[[1]]

	newRec['ncbiTaxon'] <- paste(tax[1:2], collapse = ' ')
	newRec['ncbiTaxon_ssp'] <- tax[3]
	newRec['ncbiTaxonID'] <- taxonTable[which(taxonTable$accession == addedAccessions[i, 'accession']), 'taxonID']
	newRec['repdbTaxon'] <- addedAccessions[i, 'repdbTaxon']
	newRec['locus'] <- addedAccessions[i, 'locus']
	newRec['accession'] <- addedAccessions[i, 'accession']
	newRec['filename'] <- 'downloaded, added'
	newRec['seqLength'] <- addedAccessions[i, 'seqLength']
	newRec['keep'] <- TRUE
	
	dat_repDB[[spInd]] <- rbind.data.frame(dat_repDB[[spInd]], newRec)

}




# add higher taxonomy
require(taxize)

classif <- data.table::fread('~/Dropbox/Oz_Crown_Ages/species_family.csv', data.table = FALSE) # Pyron's taxonomy
classif$Genus <- sapply(classif$Species, function(x) strsplit(x, '_')[[1]][1])
classif$Species <- gsub('_', ' ', classif$Species)

missingTaxa <- c()
for (i in 1:length(dat_repDB)) {
	
	genus <- strsplit(dat_repDB[[i]][1, 'repdbTaxon'], '\\s+')[[1]][1]
	dat_repDB[[i]]$genus <- genus
	
	if (dat_repDB[[i]][1, 'repdbTaxon'] %in% classif$Species) {
		ind <- which(classif$Species == dat_repDB[[i]][1, 'repdbTaxon'])
		dat_repDB[[i]]$family <- classif[ind, 'Family']
	} else if (genus %in% classif$Genus) {
		ind <- which(classif$Genus == genus)[1]
		dat_repDB[[i]]$family <- classif[ind, 'Family']
	} else {
		missingTaxa <- c(missingTaxa, names(dat_repDB)[i])
	}
}

missingFamily <- character(length(missingTaxa))
for (i in 1:length(missingTaxa)) {
	
	cat(i, '\n')
	Sys.sleep(time = 0.25)
	ncbi <- get_uid(missingTaxa[i], ask=F, messages=F)
	if (attributes(ncbi)$match == 'found') {
		zz <- classification(ncbi)
		missingFamily[i] <- zz[[1]][which(zz[[1]]$rank == 'family'), 'name']
	}
}

missingTaxa[which(missingFamily == '')]
missingFamily[which(missingTaxa == 'Fowlea asperrimus')] <- 'Colubridae'
missingFamily[which(missingTaxa == 'Fowlea punctulatus')] <- 'Colubridae'
missingFamily[which(missingTaxa == 'Fowlea schnurrenbergeri')] <- 'Colubridae'
missingFamily[which(missingTaxa == 'Kladirostratus acutus')] <- 'Psammophiidae'
missingFamily[which(missingTaxa == 'Pseudagkistrodon rudis')] <- 'Colubridae'
missingFamily[which(missingTaxa == 'Tropicagama temporalis')] <- 'Agamidae'

# add in higher taxonomy for taxa that had to be manually resolved.
for (i in 1:length(missingTaxa)) {
	spInd <- which(names(dat_repDB) == missingTaxa[i])
	dat_repDB[[spInd]]$family <- missingFamily[i]
}

# confirm that all elements have the same columns
table(sapply(dat_repDB, ncol))

# what's the distribution of number of loci per species?
table(sapply(dat_repDB, nrow))
table(sapply(dat_repDB, nrow) > 1)
table(sapply(dat_repDB, nrow) > 10)
table(sapply(dat_repDB, nrow) > 20)

# Now reorganize by locus and write fasta files
dat_locus <- do.call(rbind, dat_repDB)
rownames(dat_locus) <- NULL

# deal with a few last modifications
# change all Chionactis occipitalis to Sonora occipitalis
# change all Parias flavomaculatus to Trimeresurus flavomaculatus
# change taxon for RAG1 accession JN112619 to Anolis frenatus
# change taxon for RAG1 accession EU402848 to Lampropeltia californiae
# change taxon for 16S accession AF420758 to Ptychoglossus brevifrontalis
# change taxon for 16S accession DQ990972 to Xenosaurus platyceps

# in adding new accessions, should check that that taxon is not already somewhere with another taxon name

# xx <- 'Potamites juruazensis'
# matchToReptileDB(xx, considerSynonyms = TRUE)
# matchToReptileDB(xx, considerSynonyms = TRUE, interactive=T)

dat_locus <- split(dat_locus, dat_locus$locus)
sort(sapply(dat_locus, nrow))

table(sapply(dat_locus, function(x) anyDuplicated(x$accession)) > 0)
for (i in 1:length(dat_locus)) {
	if (anyDuplicated(dat_locus[[i]]$accession)) {
		message('\t', names(dat_locus)[i])
		dups <- dat_locus[[i]][which(duplicated(dat_locus[[i]]$accession) == T), 'accession']
		for (j in 1:length(dups)) {
			message('\t\t', dups[j], ':')
			message('\t\t\t', paste(which(dat_locus[[i]]$accession == dups[j]), collapse = ' '))			
		}
	}
}


fastaDir <- '~/Dropbox/Oz_Crown_Ages/fasta_March2020'



for (i in 1:length(dat_locus)) {
	
	message('\t', i, '\t', names(dat_locus)[i])
	
	locusList <- vector('list', length = nrow(dat_locus[[i]]))
	names(locusList) <- dat_locus[[i]]$accession
	
	# read in sequence data
	uniqueFiles <- unique(dat_locus[[i]]$filename)
	uniqueFiles <- uniqueFiles[!grepl('download', uniqueFiles)]
	uniqueFileList <- vector('list', length(uniqueFiles))
	names(uniqueFileList) <- uniqueFiles
	for (j in 1:length(uniqueFiles)) {
		
		f <- readDNAStringSet(uniqueFiles[j])
		# rename seqs
		seqNames1 <- setdiff(names(f), grep('@', names(f), value = TRUE))
		seqNames2 <- grep('@', names(f), value = TRUE)
		if (length(seqNames1) > 0) {
			seqNames1Acc <- strsplit(seqNames1, '\\s+')
			seqNames1Acc <- sapply(seqNames1Acc, function(x) x[1])
			seqNames1Acc <- gsub('\\.[0-9]$', '', seqNames1Acc)
			names(f)[sapply(seqNames1, function(x) which(names(f) == x))] <- seqNames1Acc
		}
		if (length(seqNames2) > 0) {			
			qq <- strsplit(seqNames2, '@') # first part is NCBI taxonID, second part is sample name/locus, and taxon name (maybe)
			seqNames2Acc <- unlist(lapply(qq, function(x) strsplit(x[2], '\\s+')[[1]][1]))
			names(f)[sapply(seqNames2, function(x) which(names(f) == x))] <- seqNames2Acc
		}
		
		uniqueFileList[[j]] <- f
	}
	
	# for each record, identify the appropriate sequence data file, and pull out the relevant sequence
	for (j in 1:nrow(dat_locus[[i]])) {
		rec <- dat_locus[[i]][j,]
		ind <- which(names(uniqueFileList) == rec$filename)
		if (length(ind) == 0 & !grepl('download', rec$filename)) stop()
		if (grepl('download', rec$filename)) {
			if (length(which(newSeqs[, 'accession'] == rec$accession & newSeqs[, 'locus'] == rec$locus)) > 0) {
				seq <- newSeqs[which(newSeqs[, 'accession'] == rec$accession & newSeqs[, 'locus'] == rec$locus), 'seq']
			} else if (length(which(addedAccessions[, 'accession'] == rec$accession & addedAccessions[, 'locus'] == rec$locus)) > 0) {
				seq <- addedAccessionsSeq[[which(addedAccessions[, 'accession'] == rec$accession & addedAccessions[, 'locus'] == rec$locus)]]
			} else {
				stop()
			}
		} else {
	
			if (rec$accession %in% names(uniqueFileList[[ind]])) {
				seq <- uniqueFileList[[ind]][rec$accession,]
			} else if (paste0(rec$accession, '_', rec$locus) %in% names(uniqueFileList[[ind]])) {
				seq <- uniqueFileList[[ind]][paste0(rec$accession, '_', rec$locus),]
			} else {
				stop()
			}
		}
		seq <- as.character(seq)
		names(seq) <- NULL	
		locusList[[j]] <- seq
	}
	names(locusList) <- dat_locus[[i]]$accession
	
	# write to file
	fn <- paste0(fastaDir, '/', names(dat_locus)[i], '.fa')
	for (j in 1:length(locusList)) {
		write(paste0('>', names(locusList)[j]), file = fn, append = ifelse(j == 1, FALSE, TRUE))
		write(locusList[[j]], file = fn, append = TRUE)
	}
	
}



## Write metadata file to the same directory
masterdat <- do.call(rbind, dat_locus)
rownames(masterdat) <- NULL

# some columns can be left out
masterdat <- masterdat[, setdiff(colnames(masterdat), 'keep')]

# sort
masterdat <- masterdat[with(masterdat, order(family, genus, repdbTaxon, locus)), ]
rownames(masterdat) <- NULL

head(masterdat)

write.csv(masterdat, '~/Dropbox/Oz_Crown_Ages/fasta_March2020/metatable.csv', row.names = FALSE)



# example
masterdat[which(masterdat$ncbiTaxon == 'Lampropeltis getula'), ]

# There are 29 (or 0.4%) of taxa that are not accepted reptileDB names. Example Corallus hortulanus

xx <- 'Corallus hortulanus'
matchToReptileDB(xx, considerSynonyms = TRUE)
matchToReptileDB(xx, considerSynonyms = TRUE, interactive=T)

masterdat[which(masterdat$repdbTaxon == 'Corallus hortulanus'), ]
