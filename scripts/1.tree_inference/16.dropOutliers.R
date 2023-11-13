# drop outliers
#	- drop samples that were embargoed
# 	- drop samples that seem out of place taxonomically
# 	- drop samples where there are multiple seqs for a single taxon
# 	- drop samples that are on the end of suspiciously long branches
# 	- drop samples that do not have any locus overlap with the rest of their genus
# 	- drop samples that have been flagged as rogues according to quartet sampling (cutoff of 0.02)


outfile <- 'samplesToDrop_all.csv'

require(Biostrings)

# dataset-wide samples to remove:
## embargoed:
embargoed <- c('K_KHOSROW', 'TJC1639', 'TJC1641', 'MBUR394', 'No_spec_No_KC_40', 'SVN_0727')

rogueFile <- '~/squam2020/roguesFromGeneTrees.csv'
overlapFile <- '~/squam2020/locusOverlap_toDrop.csv'
outlierFile <- '~/squam2020/samplesToDrop_new.csv'
flaggingFile <- '~/squam2020/bestGeneTrees/Squamata_clusterTreeReport.csv'

alnDir <- '~/squam2020/finalAlignments_June2020'

outfile <- 'samplesToDrop_all.csv'


rogues <- read.csv(rogueFile, stringsAsFactors = FALSE)
overlapOutliers <- read.csv(overlapFile, stringsAsFactors = FALSE)[,1]
outliers <- read.csv(outlierFile, stringsAsFactors = FALSE)
clusterReport <- read.csv(flaggingFile, stringsAsFactors = FALSE)
dropAsDup <- clusterReport[which(clusterReport$dropAsDup == TRUE), ]

# formatting issues
outliers$taxa <- gsub('\\s+', '', outliers$taxa)
outliers$cluster <- gsub('\\s+', '', outliers$cluster)
grep('\\s', outliers$taxa, value=T)
grep('\\−', outliers$taxa, value=T)
outliers$taxa <- gsub('\\−', '-', outliers$taxa)
grep('\\−', outliers$taxa, value=T)

grep('\\−', rogues$qf, value=T)

# were any records marked to be kept rather than discarded?
unique(outliers$reason)
keepDup <- outliers[which(outliers$reason == 'keepDup'),]

setwd(alnDir)

# read in table that allows us to interpret tip labels
metaTableFile <- '~/Dropbox/Oz_Crown_Ages/fasta_March2020/metatable2.csv'
if (file.exists(metaTableFile)) {
	taxonTable <- read.csv(metaTableFile, stringsAsFactors=FALSE)
} else {
	metaTableFile <- '~/squam2020/metatable2.csv'
	taxonTable <- read.csv(metaTableFile, stringsAsFactors=FALSE)
}
head(taxonTable)

clusterTag <- '_onedirection.trimmedFiltered_LINS1_trimmed.aln|_onedirection.filtered_LINS1_trimmed.aln'

alnFiles <- list.files(pattern = '\\.aln$')

outlierList <- vector('list', length(alnFiles))
names(outlierList) <- gsub(clusterTag, '', alnFiles)

nSamples <- numeric(length(alnFiles))
names(nSamples) <- gsub(clusterTag, '', alnFiles)

for (i in 1:length(alnFiles)) {
	
	gene <- gsub(clusterTag, '', alnFiles[i])
	message('\t', gene)
	
	aln <- readDNAStringSet(alnFiles[i])
	nSamples[i] <- length(aln)
	if (!all(names(aln) %in% taxonTable[which(taxonTable$locus == gene), 'accession'])) stop('not all accessions found!')
	gensp <- taxonTable[which(taxonTable$locus == gene),]
	gensp <- sapply(names(aln), function(x) gensp[which(gensp$accession == x), 'repdbTaxon'])
	
	# if there were any duplicates, pull those out to be flagged
	locusDropAsDup <- dropAsDup[which(dropAsDup$locus == gene), 'accession']
	if (gene %in% keepDup$cluster & length(locusDropAsDup) > 0) {
		locusDropAsDup <- setdiff(locusDropAsDup, keepDup[which(keepDup == gene), 'taxa'])
	}
		
	# embargoed
	locusEmbargoed <- unlist(sapply(embargoed, function(x) grep(x, names(aln), value = TRUE)), use.names = FALSE)
	
	# rogues
	locusRogues <- rogues[which(rogues$gene == gene), 'qf']
	
	# genus overlap
	locusOverlaps <- intersect(overlapOutliers, gensp)
	if (length(locusOverlaps) > 0) {
		locusOverlaps <- names(gensp)[unlist(lapply(locusOverlaps, function(x) which(gensp == x)))]
	}
	
	# taxonomic
	## locus specific
	locusOutliers <- outliers[which(outliers$cluster == gene), 'taxa']
	# naming may not always be perfect, double-check
	if (length(locusOutliers) > 0) {
		locusOutliers <- grep(paste0(locusOutliers, collapse = '|'), names(aln), value = TRUE)
	}
	
	## dataset-wide
	anyOutliers <- outliers[which(outliers$cluster == 'any'), 'taxa']
	if (length(anyOutliers) > 0) {
		anyOutliers <- grep(paste0(anyOutliers, collapse = '|'), names(aln), value = TRUE)
	}
	
	allOutliers <- c(locusEmbargoed, locusRogues, locusOverlaps, locusOutliers, anyOutliers, locusDropAsDup)
	allOutliers <- sort(unique(allOutliers))
	outlierMat <- matrix(0, nrow = length(allOutliers), ncol = 7)
	colnames(outlierMat) <-	c('locus', 'sample', 'embargoed', 'rogues', 'overlapOutliers', 'treeOutliers', 'duplicates')
#	outlierMat <- as.data.frame(outlierMat, stringsAsFactors = FALSE)
	outlierMat[, 1] <- gene
	outlierMat[, 2] <- allOutliers
	outlierMat[, 3] <- as.numeric(allOutliers %in% locusEmbargoed)
	outlierMat[, 4] <- as.numeric(allOutliers %in% locusRogues)
	outlierMat[, 5] <- as.numeric(allOutliers %in% locusOverlaps)
	outlierMat[, 6] <- as.numeric(allOutliers %in% c(locusOutliers, anyOutliers))
	outlierMat[, 7] <- as.numeric(allOutliers %in% locusDropAsDup)

	 outlierList[[i]] <- outlierMat	

}

sort(sapply(outlierList, nrow))
sort(sapply(outlierList, nrow) / nSamples)
range(sapply(outlierList, nrow) / nSamples)



# For each locus, remove the flagged samples and write new fasta alignments

# these will be written to the same folder as the current alignments

for (i in 1:length(alnFiles)) {
	
	gene <- gsub(clusterTag, '', alnFiles[i])
	message('\t', gene)
	
	aln <- readDNAStringSet(alnFiles[i])

	if (nrow(outlierList[[i]]) > 0) {
		
		if (length(intersect(names(aln), outlierList[[i]][, 'sample'])) != length(outlierList[[i]][, 'sample'])) stop('not perfect overlap.')
		aln2 <- aln[setdiff(names(aln), outlierList[[i]][, 'sample']),]
	} else {
		aln2 <- aln
	}
	
	fn <- gsub('\\.aln', '_outliersRemoved.aln', alnFiles[i])
	
	for (j in 1:length(aln2)) {
		write(paste0('>', names(aln2)[j]), file = fn, append = ifelse(j == 1, FALSE, TRUE))
		write(as.character(aln2[[j]]), file = fn, append = TRUE)
	}
}







