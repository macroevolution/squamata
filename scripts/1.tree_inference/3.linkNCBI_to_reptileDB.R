# Match up NCBI taxonomy with ReptileDB

# Identify accessions where taxonomy can be mapped 1-to-1 to reptileDB
# Identify accessions shared with Tonini and Zaher that need some manual updating
# Identify NCBI taxon names that are not matching to anything and set those aside for manual checking


require(pbapply)

genomTable <- read.csv('~/Dropbox/Oz_Crown_Ages/squamate_phylogenomics_v11_ncbi.csv', stringsAsFactors=FALSE)

repDBdir <- '~/Dropbox/Oz_Crown_Ages/reptileDB/'

taxonTableFile <- '~/Dropbox/Oz_Crown_Ages/Lepidosauria_8504.table'
taxonTable <- data.table::fread(taxonTableFile, data.table=FALSE)
colnames(taxonTable) <- c('V1', 'taxonID', 'locus', 'accession', 'taxon', 'definition', 'title')
head(taxonTable)

# test. Does the table correspond to the clusters?
clusterDir <- '~/Dropbox/Oz_Crown_Ages/Squamata_baitClusters_percIdentity50'

# put together list of accessions, taxa, descriptions per cluster
setwd(clusterDir)
clusters <- list.files(pattern = '\\.fa$')
clusterList <- vector('list', length(clusters))
names(clusterList) <- gsub('\\.fa', '', clusters)
for (i in 1:length(clusters)) {
	message('\t', clusters[i])
	seq <- ape::read.FASTA(clusters[i])
	message('\t\t', 100*length(intersect(taxonTable$accession, names(seq)))/length(seq), '% in table.\n')
}

# yes! Good!
rm(clusterList)

# Here, we will pull together all accessions in pyphlawd clusters that we are employing.
# For genes that appeared in multiple clusters, we will read in the separated gene fasta files rather than the original cluster fasta files 
clusterDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/fastaClusters_withNonNCBI'
clusterDir <- '~/Dropbox/Oz_Crown_Ages/fastaClusters_withNonNCBI'

# put together list of accessions, taxa, descriptions per cluster
setwd(clusterDir)
clusters <- list.files(pattern = '\\.fa$')
clusters <- sort(grep('cluster', clusters, value = TRUE))
loci <- gsub('(.+)_(cluster\\d+)', '\\1', clusters)
loci <- gsub('(cluster\\d+)_(.+)', '\\2', loci)
loci <- gsub('\\.fa$', '', loci)
loci <- (sort(unique(loci)))

clusterNums <- gsub('(.+)_(cluster\\d+)', '\\2', clusters)
clusterNums <- gsub('(cluster\\d+)_(.+)', '\\1', clusterNums)
clusterNums <- gsub('\\.fa$', '', clusterNums)

clusterDups <- unique(clusterNums[which(duplicated(clusterNums) == TRUE)])

clusterFiles <- clusters[!sapply(clusterNums, function(x) x %in% clusterDups)]

for (i in 1:length(clusterDups)) {

	grep(clusterDups[i], clusters, value = TRUE)
	f <- list.files(pattern = clusterDups[i])
	f <- grep('\\.fa\\.', f, value = TRUE)
	if (length(f) > 0) {
		f <- f[sapply(gsub('(.+)_(cluster\\d+)(\\.fa)\\.(.+)', '\\4', f), function(x) x %in% loci)]
		clusterFiles <- c(clusterFiles, f)
	} else {
		clusterFiles <- c(clusterFiles, grep(clusterDups[i], clusters, value = TRUE))
	}
}

clusterFiles <- sort(unique(clusterFiles))


clusterList <- vector('list', length(clusterFiles))
names(clusterList) <- clusterFiles


for (i in 1:length(clusterFiles)) {
	message('\t', clusterFiles[i])
	seq <- ape::read.FASTA(clusterFiles[i])
	
	if (!is.null(seq)) {
	
		locus <- gsub('_cluster[0-9]+\\.fa(\\..+)?$', '', clusterFiles[i])
		locus <- gsub('^cluster[0-9]+_|\\.fa$', '', locus)

		# pull info for seqs from genbank
		seqNames <- names(seq)
		seqNames <- gsub('_COI$|_CYTB$|_ND1$|_ND2$|_ND4$|_rRNA_12S$|_rRNA_16S$', '', seqNames)
		tmp <- taxonTable[taxonTable$accession %in% seqNames, c('taxon', 'taxonID', 'accession', 'definition')]
		
		# add info for seqs that were added in
		notInTable <- seqNames[which(!seqNames %in% taxonTable$accession)]
		if (length(notInTable) > 0) {
			
			# if there is a @ symbol, then taxon name precedes it. If there is no @ symbol but is sqcl, look it up in table
			for (j in 1:length(notInTable)) {
				
				if (grepl('@', notInTable[j])) {
					qq <- strsplit(notInTable[j], '@')[[1]] # first part is non-NCBI taxon name, second part is sample name, locus
					tmp <- rbind(tmp, cbind(taxon = gsub('_', ' ', qq[1]), taxonID = NA, accession = qq[2], definition = qq[2]))
				} else {
					if (gsub('SqCL_', '', notInTable[j]) %in% genomTable$sample) {
						zz <- unique(genomTable[which(genomTable$sample == gsub('SqCL_', '', notInTable[j])), 'lineage'])
						if (length(zz) > 1) stop('more than 1 match')
						acc <- notInTable[j]
						if (!grepl(paste0('_', locus, '$'), acc)) {
							acc <- paste0(acc, '_', locus)
						}
						tmp <- rbind(tmp, cbind(taxon = gsub('_', ' ', zz), taxonID = NA, accession = acc, definition = NA))
						
					} else {
						stop()
					}
				}
			}
		}
		rownames(tmp) <- 1:nrow(tmp)
		
		# add in locus name
		tmp <- cbind(tmp, locus = locus)

		# add in file name
		tmp <- cbind(tmp, filename = paste0('~/Dropbox/Oz_Crown_Ages/fastaClusters_withNonNCBI/', clusterFiles[i]))

	
		for (j in 1:ncol(tmp)) {
			tmp[,j] <- as.character(tmp[,j])
		}
	
		clusterList[[i]] <- tmp
	}
}

head(lapply(clusterList, head))


## CAND1 vs ND1: There has been mismanagement in terms of regular expression and these two have been confused. 
## CAND1 cluster 11871 contains all accessions in ND1 cluster 11871. In other words, ND1 cluster 11871 is completely redundant and can be removed. Why are some accessions in CAND1 cluster but not in ND1 cluster? Probably because those that are in both are those that contain the term CAND1, whereas the seqs that are in the CAND1 cluster only contain the term AHE-L288.

clusterList[['ND1_cluster11871.fa']] <- NULL

# this list is currently structured around the clusters. Instead structure aroung the loci.
qq <- do.call(rbind, clusterList)
rownames(qq) <- 1:nrow(qq)
qq <- split(qq, qq$locus)

head(lapply(qq, head))

clusterList <- qq
rm(qq)


# What about taxa that for some reason didn't make it into the pyphlawd clusters but that we are adding in?
missingDir <- '~/Dropbox/Oz_Crown_Ages/missingSeqs'
f <- list.files(missingDir, pattern = '\\.fa(sta)?$')
f <- setdiff(f, grep('onedirection|trimmed|lowmatch', f, value = TRUE))
# f <- list.files(missingDir, pattern = '_onedirection.trimmed.fa(sta)?$')
# f <- f[!f %in% c('RAG1_onedirection.trimmed.fa', 'rRNA_12S_onedirection.trimmed.fa', 'rRNA_16S_onedirection.trimmed.fa')]
# f <- sort(c(f, list.files(missingDir, pattern = 'rRNA_[0-9]{2}S_onedirection.fa$')))
f <- paste0(missingDir, '/', f)
missingList <- vector('list', length(f))
for (i in 1:length(f)) {
	message('\t', basename(f[i]))
	seq <- ape::read.FASTA(f[i])
	
	locus <- gsub('missing_|\\.fa$', '', basename(f[i]))
	locus <- gsub('_onedirection.trimmed', '', locus)

	if (length(seq) > 0) {
		
		# pull info for seqs from genbank
		seqNames1 <- setdiff(names(seq), grep('@', names(seq), value = TRUE))
		seqNames2 <- grep('@', names(seq), value = TRUE)
		seqNames1Acc <- strsplit(seqNames1, '\\s+')
		seqNames1Acc <- sapply(seqNames1Acc, function(x) x[1])
		seqNames1Acc <- gsub('\\.[0-9]$', '', seqNames1Acc)
		tmp <- taxonTable[taxonTable$accession %in% seqNames1Acc, c('taxon', 'taxonID', 'accession', 'definition')]
			
		# add info for seqs that were added in
		if (length(seqNames2) > 0) {
			
			qq <- strsplit(seqNames2, '@') # first part is NCBI taxonID, second part is sample name/locus, and taxon name (maybe)
			samp <- unlist(lapply(qq, function(x) gsub('SqCL_', '', x[2])))
			samp <- gsub('rRNA_12S|rRNA_16S', 'rRNA', samp)
			samp <- sapply(strsplit(samp, '_'), function(x) paste(x[- length(x)], collapse='_'))
			tax <- sapply(samp, function(x) genomTable[which(genomTable$sampl == x), 'lineage'])
			acc <- unlist(lapply(qq, function(x) strsplit(x[2], '\\s+')[[1]][1]))
			tmp <- rbind(tmp, cbind(taxon = gsub('_', ' ', tax), taxonID = NA, accession = acc, definition = unlist(lapply(qq, function(x) x[2]))))
			
			# # qq <- strsplit(seqNames2, '@') # first part is NCBI taxonID, second part is sample name/locus, and taxon name
			# qq <- sapply(qq, function(x) c(x[1], strsplit(x[2], '\\s+')[[1]]), simplify = FALSE)
			# qq <- lapply(qq, function(x) c(x[1], x[2], paste(x[3:length(x)], collapse = ' ')))
			# tmp <- rbind(tmp, cbind(taxon = unlist(lapply(qq, function(x) x[3])), taxonID = unlist(lapply(qq, function(x) x[1])), accession = unlist(lapply(qq, function(x) x[2])), definition = unlist(lapply(qq, function(x) x[2]))))
		}
		rownames(tmp) <- 1:nrow(tmp)
		
	#	if (any(!grepl(paste0('_', locus, '$'), tmp$accession) & grepl('sqcl', tmp$accession, ignore.case=TRUE))) stop()
		
		# add in locus name
		tmp <- cbind(tmp, locus = locus)
		
		# add in filename
		tmp <- cbind(tmp, filename = f[i])
		
		for (j in 1:ncol(tmp)) {
			tmp[,j] <- as.character(tmp[,j])
		}
		
		missingList[[i]] <- tmp
	} else {
		message('\t\tempty.')
	}
}

head(lapply(missingList, head))

# Add missing list to the clusterList. This will make it easier to know which loci each taxon has.
for (i in 1:length(clusterList)) {
	locus <- clusterList[[i]][1, 'locus']
	if (locus %in% unlist(lapply(missingList, function(x) x[1, 'locus']))) {
		ind <- which(unlist(lapply(missingList, function(x) x[1, 'locus'])) == locus)
		clusterList[[i]] <- rbind.data.frame(clusterList[[i]], missingList[[ind]], stringsAsFactors = FALSE)
	}
}

# standardize taxon naming: replace underscores with spaces and remove trailing spaces
for (i in 1:length(clusterList)) {
	tmp <- clusterList[[i]]$taxon
	tmp <- gsub('_', ' ', tmp)
	tmp <- gsub('\\s+', ' ', tmp)
	tmp <- gsub('(^\\s)|(\\s$)', '', tmp)
	clusterList[[i]]$taxon <- tmp
}

########################################
# There is a mismanagement of accessions where accessions that are CAND1 may be listed under ND1 and vice versa. We will identify those problems and rectify here.

## THIS IS NOW HANDLED ABOVE BY SIMPLY DROPPING THE ND1 CLUSTER THAT WAS IN FACT 100% REDUNDANT WITH THE CAND1 CLUSTER
## NO RECORDS SHOULD BE APPEARING AS MODIFIED.

cand1Regex <- '(cullin-associated and neddylation-dissociated protein 1)|(CAND1)'
nd1Regex <- '((?<!CA)ND1)|(dehydrogenase subunit 1)|(mitochondrion complete genome)|(mitochondrial DNA complete genome)|(mitochondrion partial genome)'


nd1Clusters <- grep('^nd1', names(clusterList), ignore.case = TRUE, value = TRUE)
cand1Clusters <- grep('cand1', names(clusterList), ignore.case = TRUE, value = TRUE)

for (i in 1:length(nd1Clusters)) {
	
	for (j in 1:nrow(clusterList[[nd1Clusters[i]]])) {
		
		table(clusterList[[nd1Clusters[i]]][,'locus'])
		
		tmp <- clusterList[[nd1Clusters[i]]][j,]

		# If a SqCL record, or if the seq contained multiple gene regions and was split, identifying info will be in accession, not definition
		if (grepl('CAND1|AHE-L288|(?<!CA)ND1', tmp$accession, ignore.case = TRUE, perl = TRUE)) {
			cand1Found <- grepl('CAND1|AHE-L288', tmp$accession, ignore.case = TRUE)
			nd1Found <- grepl('(?<!CA)ND1', tmp$accession, ignore.case = TRUE, perl = TRUE)				
		} else {
			cand1Found <- grepl(cand1Regex, tmp$definition, ignore.case = TRUE)
			nd1Found <- grepl(nd1Regex, tmp$definition, ignore.case = TRUE, perl = TRUE)
		}
		
		if (cand1Found & clusterList[[nd1Clusters[i]]][j, 'locus'] == 'ND1') {
			clusterList[[nd1Clusters[i]]][j, 'locus'] <- 'CAND1'
			message('\t\tfixed.')
		}
		if (nd1Found & clusterList[[nd1Clusters[i]]][j, 'locus'] == 'CAND1') {
			clusterList[[nd1Clusters[i]]][j, 'locus'] <- 'ND1'
			message('\t\tfixed.')
		}
	}
}

for (i in 1:length(cand1Clusters)) {
	
	for (j in 1:nrow(clusterList[[cand1Clusters[i]]])) {
		
		table(clusterList[[cand1Clusters[i]]][,'locus'])
		
		tmp <- clusterList[[cand1Clusters[i]]][j,]

		# If a SqCL record, or if the seq contained multiple gene regions and was split, identifying info will be in accession, not definition
		if (grepl('CAND1|AHE-L288|(?<!CA)ND1', tmp$accession, ignore.case = TRUE, perl = TRUE)) {
			cand1Found <- grepl('CAND1|AHE-L288', tmp$accession, ignore.case = TRUE)
			nd1Found <- grepl('(?<!CA)ND1', tmp$accession, ignore.case = TRUE, perl = TRUE)				
		} else {
			cand1Found <- grepl(cand1Regex, tmp$definition, ignore.case = TRUE)
			nd1Found <- grepl(nd1Regex, tmp$definition, ignore.case = TRUE, perl = TRUE)
		}
		
		if (cand1Found & clusterList[[cand1Clusters[i]]][j, 'locus'] == 'ND1') {
			clusterList[[cand1Clusters[i]]][j, 'locus'] <- 'CAND1'
			message('\t\tfixed.')
		}
		if (nd1Found & clusterList[[cand1Clusters[i]]][j, 'locus'] == 'CAND1') {
			clusterList[[cand1Clusters[i]]][j, 'locus'] <- 'ND1'
			message('\t\tfixed.')
		}
	}
}




######
# saveRDS(clusterList, '~/Dropbox/Oz_Crown_Ages/clusterList_4Feb2020.rds')
# Load clusterList from saved file
# clusterList <- readRDS('~/Dropbox/Oz_Crown_Ages/clusterList_4Feb2020.rds')
######

# what are the unique NCBI taxa?
# ncbiTaxa <- sort(unique(taxonTable$taxon))
ncbiTaxa <- unlist(lapply(clusterList, function(x) x$taxon))
ncbiTaxa <- gsub('_', ' ', ncbiTaxa)
ncbiTaxa <- gsub('\\s+', ' ', ncbiTaxa)
ncbiTaxa <- gsub('(^\\s)|(\\s$)', '', ncbiTaxa)
ncbiTaxa <- sort(unique(ncbiTaxa))
# ncbiTaxa <- union(ncbiTaxa, taxonTable$taxon)
length(ncbiTaxa)

# Read in Reptile Database data
acceptedNamesTable <- read.csv(paste0(repDBdir, 'acceptedTaxonTable.csv'), stringsAsFactors = FALSE)
subspeciesTable <- read.csv(paste0(repDBdir, 'subspeciesTable.csv'), stringsAsFactors = FALSE)
synonymTable <- read.csv(paste0(repDBdir, 'synonymTable.csv'), stringsAsFactors = FALSE)

# limit synonyms to post-1980
synonymTable <- synonymTable[which(synonymTable$year >= 1980 | is.na(synonymTable$year)),]

source('~/Dropbox/Oz_Crown_Ages/reptileDB/matchToReptileDB.R')


# Go through all NCBI taxon names, and associate a ReptileDB name only when there is a direct mapping between the two.
# We will:
# 	Look for exact matching with accepted names and subspecies of accepted names
# 	Look at synonyms that postdate 1980, even if there is already an exact match with an accepted name.
# 	Place in results table: 
# 							NCBI taxon name and ID
# 							associated repDB taxon name if it was found
# 							accession and name of locus
# 							reason: exact mapping, multiple matches, etc

resList <- vector('list', length(ncbiTaxa))
names(resList) <- ncbiTaxa

for (i in 1:length(ncbiTaxa)) {
	
	message('\t', i, '\t', ncbiTaxa[i])

	# which loci is this taxon present in, and what are the accessions?
	loci <- which(sapply(clusterList, function(x) ncbiTaxa[i] %in% x$taxon) == TRUE)
	loci <- lapply(clusterList, function(x) x[which(x$taxon == ncbiTaxa[i]),])
	loci <- loci[which(sapply(loci, nrow) > 0)]
	loci <- do.call(rbind, loci)
	
	match <- matchToReptileDB(ncbiTaxa[i], considerSynonyms = TRUE, resolveToSubSpecies = FALSE)
	
	resTable <- matrix(nrow = nrow(loci), ncol = 9)
	colnames(resTable) <- c('ncbiTaxon', 'ncbiTaxon_ssp', 'ncbiTaxonID', 'repdbTaxon', 'repdbTaxon_ssp', 'locus', 'accession', 'reason', 'filename')
	
	# fill in NCBI taxon info
	resTable[, 'ncbiTaxon'] <- match['ncbiTaxon']
	resTable[, 'ncbiTaxon_ssp'] <- match['ncbiTaxon_ssp']
	resTable[, 'ncbiTaxonID'] <- taxonTable[which(taxonTable$taxon == ncbiTaxa[i])[1], 'taxonID']
	resTable[, 'locus'] <- loci$locus
	resTable[, 'accession'] <- loci$accession
	resTable[, 'repdbTaxon'] <- match['repdbTaxon']
	resTable[, 'repdbTaxon_ssp'] <- match['repdbTaxon_ssp']
	resTable[, 'reason'] <- match['reason']
	resTable[, 'filename'] <- loci$filename
	
	# occasionally, there may be duplicate entries (such as due to corrected error of CAND1 vs ND1)
	resTable <- unique(resTable)
	
	for (j in 1:ncol(resTable)) {
		resTable[,j] <- as.character(resTable[,j])
	}

	resList[[i]] <- resTable
	
}


byNCBItaxon <- as.data.frame(do.call(rbind, resList), stringsAsFactors=FALSE)
head(byNCBItaxon)

# how many don't have a reptileDB taxon name?
singleRep <- lapply(split(byNCBItaxon, apply(byNCBItaxon, 1, function(x) paste(x[1:2], collapse = ' '))), function(x) x[1,])
table(is.na(sapply(singleRep, function(x) x$repdbTaxon)))

# what is the breakdown of reasons?
table(sapply(singleRep, function(x) x$reason), useNA = 'always') # for taxa
table(byNCBItaxon$reason, useNA = 'always')



###################
## Tonini and Zaher


# load accession numbers by species from Tonini 2016 and Zaher 2019.
## We will use these to update our own effort at taxonomic reconciliation.
toniniTable <- read.csv('~/Dropbox/Oz_Crown_Ages/Tonini2016/squam_shl_new_GenBank.csv', stringsAsFactors=FALSE)
zaherTable <- gdata::read.xls('~/Dropbox/Oz_Crown_Ages/Zaher2019/ZaherTable_S3.xlsx', stringsAsFactors=FALSE)


toniniTable <- toniniTable[which(toniniTable$Genes > 0 & toniniTable$Taxon %in% c('Amphisbaenia', 'Rhynchocephalia', 'Sauria', 'Serpentes')),]
colnames(toniniTable)[grep('X12S', colnames(toniniTable))] <- 'rRNA_12S'
colnames(toniniTable)[grep('X16S', colnames(toniniTable))] <- 'rRNA_16S'
colnames(toniniTable)[grep('CMOS', colnames(toniniTable))] <- 'cmos'
colnames(toniniTable)[grep('X12S', colnames(toniniTable))] <- 'rRNA_12S'
dim(toniniTable)
tail(toniniTable)

zaherTable <- zaherTable[which(zaherTable$family != ''),]
zaherTable$terminal <- gsub('_', ' ', zaherTable$terminal)
head(zaherTable)


byNCBItaxon$toniniTaxon <- NA
byNCBItaxon$zaherTaxon <- NA

for (i in 1:nrow(toniniTable)) {
	
	acc <- toniniTable[i, 4:ncol(toniniTable)]
	acc <- setdiff(as.character(acc), '')
	toniniTaxon <- toniniTable[i, 'Species']
	
	if (any(acc %in% byNCBItaxon$accession)) {
				
		# for each of the accessions that are shared, how does the taxon name in Tonini compare?
		byNCBItaxon[byNCBItaxon$accession %in% acc, 'toniniTaxon'] <- toniniTaxon
		
	}	
}

# where we did fill in a taxon name from Tonini, do we find:
# 		difference in resolved taxon name?
# 		Tonini providing a name when we failed to find a match?

# How many accessions are shared?
table(is.na(byNCBItaxon$toniniTaxon))

# 	difference in resolved taxon name?
tmp <- byNCBItaxon[which(!is.na(byNCBItaxon$toniniTaxon) & !is.na(byNCBItaxon$repdbTaxon)),]
table(tmp$toniniTaxon != tmp$repdbTaxon)

tmp <- tmp[which((tmp$toniniTaxon != tmp$repdbTaxon) == TRUE),]
length(split(tmp, tmp$ncbiTaxon)) # number of unique ncbi taxa
head(tmp)

# 	Tonini providing a name when we failed to find a match?
tmp <- byNCBItaxon[which(!is.na(byNCBItaxon$toniniTaxon) & is.na(byNCBItaxon$repdbTaxon)),]
nrow(tmp) # number of accessions
length(split(tmp, tmp$ncbiTaxon)) # number of unique ncbi taxa
head(tmp, 30)


## Do the same with Zaher

for (i in 1:nrow(zaherTable)) {
	
	acc <- zaherTable[i, 10:ncol(zaherTable)]
	acc <- setdiff(as.character(acc), c('', '#N/A'))
	zaherTaxon <- zaherTable[i, 'terminal']
	
	if (any(acc %in% byNCBItaxon$accession)) {
	
		# for each of the accessions that are shared, how does the taxon name in Tonini compare?
		byNCBItaxon[byNCBItaxon$accession %in% acc, 'zaherTaxon'] <- zaherTaxon
	}	
}


# where we did fill in a taxon name from Zaher, do we find:
# 		difference in resolved taxon name?
# 		tonini providing a name when we failed to find a match?

# 	difference in resolved taxon name?
tmp <- byNCBItaxon[which(!is.na(byNCBItaxon$zaherTaxon) & !is.na(byNCBItaxon$repdbTaxon)),]
table(tmp$zaherTaxon != tmp$repdbTaxon)

tmp <- tmp[which((tmp$zaherTaxon != tmp$repdbTaxon) == TRUE),]
head(tmp)

# 	Zaher providing a name when we failed to find a match?
tmp <- byNCBItaxon[which(!is.na(byNCBItaxon$zaherTaxon) & is.na(byNCBItaxon$repdbTaxon)),]
nrow(tmp) # number of accessions
length(split(tmp, tmp$ncbiTaxon)) # number of unique ncbi taxa
head(tmp, 30)


# For both Tonini and Zaher, add a matched reptileDB name, in case their taxon names are now outdated. This will also help reconcile their names with my matches from ReptileDB. 
## If there is ambiguity in terms of which name to match to, consider the year associated with the names. 
## Tonini built up their taxonomy based on Reptile Database circa 2015.
## Zaher built up their taxonomy based partially on Reptile Database circa 2017.


byNCBItaxon$toniniMatch <- NA
byNCBItaxon$zaherMatch <- NA


for (i in 1:nrow(byNCBItaxon)) {

	if (!is.na(byNCBItaxon[i, 'toniniTaxon'])) {
		if (byNCBItaxon[i, 'toniniTaxon'] != byNCBItaxon[i, 'repdbTaxon'] | is.na(byNCBItaxon[i, 'repdbTaxon'])) {
			
			match <- matchToReptileDB(byNCBItaxon[i, 'toniniTaxon'], returnMultipleMatches = TRUE)
			
			if (is.data.frame(match)) {
				# if one of the matched taxa is from 2015 or later, then don't ignore the ambiguity
				if (any(match$year >= 2015)) {
					byNCBItaxon[i, 'toniniMatch'] <- NA
				} else {
					# if all predate 2015, then determination has already been made
					byNCBItaxon[i, 'toniniMatch'] <- byNCBItaxon[i, 'toniniTaxon']
				}
			} else {
				byNCBItaxon[i, 'toniniMatch'] <- match['repdbTaxon']
			}
		} else {
			byNCBItaxon[i, 'toniniMatch'] <- byNCBItaxon[i, 'repdbTaxon']
		}	
	}

	if (!is.na(byNCBItaxon[i, 'zaherTaxon'])) {
		if (byNCBItaxon[i, 'zaherTaxon'] != byNCBItaxon[i, 'repdbTaxon'] | is.na(byNCBItaxon[i, 'repdbTaxon'])) {
			
			match <- matchToReptileDB(byNCBItaxon[i, 'zaherTaxon'], returnMultipleMatches = TRUE)
			
			if (is.data.frame(match)) {
				# if one of the matched taxa is from 2015 or later, then don't ignore the ambiguity
				if (any(match$year >= 2017)) {
					byNCBItaxon[i, 'zaherMatch'] <- NA
				} else {
					# if all predate 2017, then determination has already been made
					byNCBItaxon[i, 'zaherMatch'] <- byNCBItaxon[i, 'zaherTaxon']
				}
			} else {
				byNCBItaxon[i, 'zaherMatch'] <- match['repdbTaxon']
			}
		} else {
			byNCBItaxon[i, 'zaherMatch'] <- byNCBItaxon[i, 'repdbTaxon']
		}	
	}
}


# which records had a taxon in Tonini, but failed at exact matching?
tmp <- byNCBItaxon[which(!is.na(byNCBItaxon$toniniTaxon) & is.na(byNCBItaxon$toniniMatch)),]
nrow(tmp)
nrow(unique(tmp[, c('ncbiTaxon', 'ncbiTaxon_ssp', 'repdbTaxon', 'repdbTaxon_ssp', 'reason', 'toniniTaxon', 'toniniMatch')]))
head(unique(tmp[, c('ncbiTaxon', 'ncbiTaxon_ssp', 'repdbTaxon', 'repdbTaxon_ssp', 'reason', 'toniniTaxon', 'toniniMatch')]), 20)

# which records had a taxon in Zaher, but failed at exact matching?
tmp <- byNCBItaxon[which(!is.na(byNCBItaxon$zaherTaxon) & is.na(byNCBItaxon$zaherMatch)),]
nrow(tmp)
head(unique(tmp[, c('ncbiTaxon', 'ncbiTaxon_ssp', 'repdbTaxon', 'repdbTaxon_ssp', 'reason', 'zaherTaxon', 'zaherMatch')]), 20)


## Can we seek out just the record mismatches that are due to identified errors (rather than taxonomy updates)
### Rationale: If we include taxonomy updates for the overlapping accessions, but don't for the remaining accessions in our dataset, 
### 	then wouldn't that lead to a problematic double track?

# Here, rather than flag anything that does not perfectly match our reptileDB match, or anything that maps to multiple potential taxa, we will flag cases where:
# -	RepDB accepted name is different genus than Tonini/Zaher genus (post matching, or pre-matching if Tonini/Zaher name is not 1-to-1 matched to RepDB)
# - Same as above but with NCBI taxon name if NCBI taxon is unresolved.
# - Tonini and Zaher disagree.


byNCBItaxon$toniniZaherFlag <- FALSE

# Flag those worth tackling

for (i in 1:nrow(byNCBItaxon)) {
	
	# what is the genus in our dataset?
	if (is.na(byNCBItaxon[i, 'repdbTaxon'])) {
		origTaxon <- byNCBItaxon[i, 'ncbiTaxon']
	} else {
		origTaxon <- byNCBItaxon[i, 'repdbTaxon']
	}	
	
	toniniTaxon <- NA
	zaherTaxon <- NA
	
	# What is the genus in the Tonini dataset?
	if (!is.na(byNCBItaxon[i, 'toniniTaxon'])) {

		if (is.na(byNCBItaxon[i, 'toniniMatch'])) {
			toniniTaxon <- byNCBItaxon[i, 'toniniTaxon']
		} else {
			toniniTaxon <- byNCBItaxon[i, 'toniniMatch']
		}	
	}


	# What is the genus in the Zaher dataset?
	if (!is.na(byNCBItaxon[i, 'zaherTaxon'])) {

		if (is.na(byNCBItaxon[i, 'zaherMatch'])) {
			zaherTaxon <- byNCBItaxon[i, 'zaherTaxon']
		} else {
			zaherTaxon <- byNCBItaxon[i, 'zaherMatch']
		}	
	}		

	if (!is.na(toniniTaxon) & !identical(origTaxon, toniniTaxon)) {
		byNCBItaxon[i, 'toniniZaherFlag'] <- TRUE
	}

	if (!is.na(zaherTaxon) & !identical(origTaxon, zaherTaxon)) {
		byNCBItaxon[i, 'toniniZaherFlag'] <- TRUE
	}
	
	if (!is.na(toniniTaxon) & !is.na(zaherTaxon)) {
		if (!identical(toniniTaxon, zaherTaxon)) {
			byNCBItaxon[i, 'toniniZaherFlag'] <- TRUE
		}
	}

}



table(byNCBItaxon$toniniZaherFlag)

tmp <- byNCBItaxon[which(byNCBItaxon$toniniZaherFlag == TRUE),]
head(tmp)

tmp$manualIdentification <- ''
tmp$comments <- ''
tmp <- tmp[with(tmp, order(ncbiTaxon, ncbiTaxon_ssp)),]

# write.csv(tmp, '~/Dropbox/Oz_Crown_Ages/toniniZaher_resolution.csv', row.names = FALSE)





sp <- 'Alopoglossus angulatus'
matchToReptileDB(sp, interactive=TRUE)
matchToReptileDB(sp, returnMultipleMatches=TRUE)


#########################
## Unmatched taxon names 

# which NCBI taxa aren't matching to anything?

tmp <- which(is.na(byNCBItaxon$repdbTaxon) & is.na(byNCBItaxon$reason))
tmp <- byNCBItaxon[tmp,]
tmp <- unique(tmp[, c('ncbiTaxon', 'ncbiTaxon_ssp', 'ncbiTaxonID')])
tmp$manualIdentification <- ''
tmp$comments <- ''

tmp <- tmp[with(tmp, order(ncbiTaxon, ncbiTaxon_ssp)),]

# write.csv(tmp, '~/Dropbox/Oz_Crown_Ages/ncbi_noMatchTaxa.csv', row.names = FALSE)



	
## SqCL
# For anything with SqCL, presumably Sonal already sorted out taxonomy, so we might be able to assume that those don't need taxonomic reconciliation.
sqClInd <- grep('sqcl', byNCBItaxon$accession, ignore.case = TRUE)

# how many rows are SqCL?
length(sqClInd)
length(sqClInd) / nrow(byNCBItaxon)

# how many that are SqCL don't have a resolved name?
length(which(is.na(byNCBItaxon$repdbTaxon) & grepl('sqcl', byNCBItaxon$accession, ignore.case = TRUE))) / length(sqClInd)

# how many sqCL TAXA don't have a resolved name?
tmp <- byNCBItaxon[which(is.na(byNCBItaxon$repdbTaxon) & grepl('sqcl', byNCBItaxon$accession, ignore.case = TRUE)),]
length(split(tmp, tmp$ncbiTaxon)) 


# how many unresolved taxa are not matched at all?
table(sapply(split(tmp, tmp$ncbiTaxon), function(x) x[1,'reason']), useNA = 'always')

unique(tmp[which(is.na(tmp$reason)), c('ncbiTaxon', 'ncbiTaxon_ssp', 'ncbiTaxonID')])

# Which taxa are multi-matches?
multi <- unique(tmp[which(tmp$reason == 'multiple matches'), c('ncbiTaxon', 'ncbiTaxon_ssp', 'ncbiTaxonID')])
nrow(multi)
head(multi, 50)

# for each of the unresolved taxa, is there a resolved subspecies that has similar gene coverage?
# If there was, then the fact that one is unresolved almost doesn't matter...
for (i in 1:nrow(multi)) {
	
	sp <- multi[i, 'ncbiTaxon']
	qq <- which(byNCBItaxon$ncbiTaxon == sp & !is.na(byNCBItaxon$ncbiTaxon_ssp))
	if (length(qq) > 0) {
		qq <- byNCBItaxon[qq,]
		altLoci <- qq[, 'locus']
		locusCoverage <- byNCBItaxon[which(byNCBItaxon$ncbiTaxonID == multi[i, 'ncbiTaxonID']), 'locus']
		message('\t', sp, ': ',length(intersect(altLoci, locusCoverage)), ' of ', length(locusCoverage))
	} else {
		message('\t', sp, ': no alternatives.')
	}
}

### Are there any accessions that are duplicated with the same gene?
table(duplicated(byNCBItaxon[,c('locus','accession')])) 



# Write file
# This file contains all accessions currently included in pyphlawd clusters that we are using and contains 1 seq per locus per taxon, but where taxon was defined according to NCBI taxonomy. The diffrence here is that it now contains reptileDB matches.

write.csv(byNCBItaxon, '~/Dropbox/squamatePhylo/2019/topologyScripts/allAccessionsMatchingToReptileDB.csv', row.names=FALSE)


