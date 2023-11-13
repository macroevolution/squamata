# flag potential outliers in alignments/trees

# Tests:
# 	- How many / what percent private alleles does sample have?
# 	- How long is the terminal branch length?
# 	- What is patristic distance to closest tip
#
# 	- What is median patristic distance to sample's own family and genus?
# 	- What is smallest median patristic distance to other families and genera?
# 	- If sample is a taxon that has subspecies in tree, do all subspecies form exclusive clade?

# Are all families and genera monophyletic?


setwd('~/squam2020/bestGeneTrees')

clusterDir <- getwd()

options(warn=1)


require(ape)
require(pbapply)
require(Biostrings)
require(rBLAST)

toniniTreeFile <- '~/Documents/pyphlawd/tonini5416.tre'
if (file.exists(toniniTreeFile)) {
	tonini <- read.tree(toniniTreeFile)
} else {
	toniniTreeFile <- '~/Downloads/Tonini_Squamata_Tree_Methods_Data/tonini5416.tre'
}
tonini <- read.tree(toniniTreeFile)
tonini$tip.label <- gsub('_', ' ', tonini$tip.label)

# read in table that allows us to interpret tip labels
taxonTableFile <- '~/Dropbox/Oz_Crown_Ages/fasta_March2020/metatable2.csv'
if (file.exists(taxonTableFile)) {
	taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
} else {
	taxonTableFile <- '~/squam2020/metatable2.csv'
	taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
}
head(taxonTable)

# these trees were already relabeled to NCBI taxon ID
keepclusters <- list.files(clusterDir, pattern = '\\.tre$|\\.treefile$')
keepclusters <- setdiff(keepclusters, grep('ncbi\\.tre', keepclusters, value=TRUE))
keepclusters <- grep('LINS1', keepclusters, value=TRUE)

clusterTag <- '_onedirection.trimmedFiltered_LINS1_trimmed.treefile|_onedirection.filtered_LINS1_trimmed.treefile'

# for BLAST
ungappedDir <- '~/Dropbox/Oz_Crown_Ages/fasta_March2020'
ungappedDir <- '~/squam2020/fasta_March2020_trimmedFiltered'
ungappedSeqFiles <- list.files(ungappedDir, pattern = '_onedirection.trimmedFiltered.fa$', full.names = TRUE)
ungappedSeqFiles <- sort(c(ungappedSeqFiles, list.files(ungappedDir, pattern = 'rRNA_1\\dS_onedirectionFiltered.fa$', full.names = TRUE)))


reportList <- vector('list', length = length(keepclusters))
names(reportList) <- gsub(clusterTag, '', keepclusters)

# new version now that we have switched to reptileDB taxonomy. 
# now, trees and alignments are labeled according to accession, and those can be converted to taxon names through the metadata table.
# the metadata table also contains higher taxonomy

for (i in 1:length(keepclusters)) {

	gene <- gsub(clusterTag, '', keepclusters[i])
	cat('\t', gene, '\n')

	tree <- read.tree(keepclusters[i])
	
	# superCRUNCH may have added * to the end of some labels. Remove those. 
	tree$tip.label <- gsub('\\*', '', tree$tip.label)

	# build up taxonomy table
	tax <- as.data.frame(matrix(nrow = Ntip(tree), ncol = 5))
	colnames(tax) <- c('accession', 'seqLength', 'family','genus','species')
	for (j in 1:length(tree$tip.label)) {
		ind <- which(taxonTable$locus == gene & taxonTable$accession == tree$tip.label[j])
		if (length(ind) != 1) stop()
		tax[j,] <- taxonTable[ind, c('accession', 'seqLength', 'family', 'genus', 'repdbTaxon')]
	}
	
	# because of the gene splitting done prior to alignment, there may now be some duplicate species.
	# we will track these and:
	# 	out of multiple seqs per species, if one is misbehaving and another isn't, we keep the the latter
	# 	and/or we keep the longer seq
	if (nrow(tax) != length(unique(tax$species))) {
		message('\t\t', sum(duplicated(tax$species)), ' duplicate taxa found')
	}
	tree$tip.label <- tax$species

	# list of taxa grouped by genus
	genusGroups <- split(tax$species, tax$genus)
	
	# list of taxa grouped by family
	familyGroups <- split(tax$species, tax$family)
	
	
	# subset tonini tree to as many of these taxa as possible
	toniniSub <- keep.tip(tonini, intersect(tonini$tip.label, tax$species))

	# what is the basal most taxon?
	xx <- which(toniniSub$edge[,1] == Ntip(toniniSub) + 1)
	xx <- toniniSub$edge[xx, 2]

	if (length(geiger::tips(toniniSub, xx[1])) < length(geiger::tips(toniniSub, xx[2]))) {
		outgroup <- geiger::tips(toniniSub, xx[1])
	} else {
		outgroup <- geiger::tips(toniniSub, xx[2])
	}

	if (any(grepl('Sphenodon', tax$species))) {
		outgroup <- grep('Sphenodon', tax$species, value=T)
	}
	
	# reroot
	if (length(outgroup) == 1) {
		tree <- ladderize(root(tree, outgroup))
	} else {
		tree <- ladderize(root(tree, node=getMRCA(tree, outgroup)))
	}

	# render ultrametric
	tree <- phytools::force.ultrametric(tree, method='nnls')

	# prepare report
	geneReport <- as.data.frame(matrix(nrow = nrow(tax), ncol = 15))
	colnames(geneReport) <- c('locus', 'accession', 'family', 'genus', 'species', 'inTonini', 'patristicOutlier', 'familyMonophyletic', 'genusMonophyletic', 'terminalEdge', 'trimmed', 'percentPrivateAlleles', 'blastFamilyPercent', 'blastGenusPercent', 'dropAsDup')

	geneReport$locus <- gene
	geneReport$accession <- tax$accession
	geneReport$family <- tax$family
	geneReport$genus <- tax$genus
	geneReport$species <- tax$species
	geneReport$inTonini <- tax$species %in% tonini$tip.label
	geneReport$patristicOutlier <- FALSE
	geneReport$familyMonophyletic <- TRUE
	geneReport$genusMonophyletic <- TRUE
	geneReport$terminalEdge <- NA
	geneReport$trimmed <- FALSE
	geneReport$percentPrivateAlleles <- 0
	geneReport$blastFamilyPercent <- NA
	geneReport$blastGenusPercent <- NA
	geneReport$dropAsDup <- FALSE

	#######
	# PATRISTIC DISTANCE BASED METRICS v2

	# For each taxon, calculate median patristic distance to each genus.
	# if smallest genus distance is the taxon's genus, good. Otherwise:
	# if distance to its own genus is smaller than distance to other genera in less than 95% of the cases, then outlier.
	# if monotypic genus, if closest family is not its own family, then outlier.
	
	# get pairwise patristic distance
	pat <- cophenetic(tree)
	diag(pat) <- NA

	for (j in 1:nrow(tax)) {
		
		focalSp <- as.character(tax[j, 'species'])
		
		meanFamilyDist <- sapply(familyGroups, function(x) median(pat[as.character(x), focalSp], na.rm=TRUE))
		meanGenusDist <- sapply(genusGroups, function(x) median(pat[as.character(x), focalSp], na.rm=TRUE))
		
		if (names(which.min(meanGenusDist)) != tax[j, 'genus']) {
			if (length(which(meanGenusDist[tax[j, 'genus']] < meanGenusDist[setdiff(names(meanGenusDist), tax[j, 'genus'])])) / length(meanGenusDist) < 0.95) {
				geneReport[j, 'patristicOutlier'] <- TRUE
			}
		
			# if genus was monotypic, check with family
			if (sum(tax$genus == tax[j, 'genus']) == 1 & names(which.min(meanFamilyDist)) != tax[j, 'genus']) {
				geneReport[j, 'patristicOutlier'] <- FALSE
			}
		}
	}
	
	########
	# MONOPHYLY CHECK
	
	# Are families monophyletic?
	for (j in 1:length(familyGroups)) {
		if (length(familyGroups[[j]]) > 1) {
			node <- getMRCA(tree, familyGroups[[j]])
			subphy <- extract.clade(tree, node = node)
			# tax[which(tax$genbankID %in% subphy$tip.label), 'species']
			if (Ntip(subphy) != length(familyGroups[[j]])) {
				geneReport[which(geneReport$family == names(familyGroups)[j]), 'familyMonophyletic'] <- FALSE
			}
		}
	}
	
	# Are genera monophyletic?
	for (j in 1:length(genusGroups)) {
		if (length(genusGroups[[j]]) > 1) {
			node <- getMRCA(tree, genusGroups[[j]])
			subphy <- extract.clade(tree, node = node)
			# tax[which(tax$genbankID %in% subphy$tip.label), 'species']
			if (Ntip(subphy) != length(genusGroups[[j]])) {
				geneReport[which(geneReport$genus == names(genusGroups)[j]), 'genusMonophyletic'] <- FALSE
			}
		}
	}

	########
	# TERMINAL BRANCH LENGTHS
	
	geneReport$terminalEdge <- sapply(1:length(tree$tip.label), function(x) {
		tree$edge.length[which(tree$edge[,2] == x)]	
	})

	# Using PyPhlawd trim_tips.py, flag tips that would be removed with a relative cutoff of 0.25, and given a 10x difference between sister lineages.
	# We will not use an absolute length filter (hence 50)
	args <- c('~/PyPHLAWD/src/trim_tips.py', keepclusters[i], 0.25, 50)

	trim <- system2('python3', args, stdout = TRUE, stderr = TRUE)
	
	trim <- grep('removing', trim, value=TRUE)
	if (length(trim) > 0) {
		for (j in 1:length(trim)) {
			if (!grepl('^removing', trim[j])) {
				trim[j] <- strsplit(trim[j], ';')[[1]][2]
			}
		}
		trim <- strsplit(trim, ' ')
		trim <- sapply(trim, function(x) x[2])
		geneReport$trimmed[geneReport$accession %in% trim] <- TRUE
	}
	
	# also use PyPhlawd trim_internal_edges.py with absolute cutoff of 1.5
	# this will catch the occasional very long branch subtending more than one tip
	args <- c('~/PyPHLAWD/src/trim_internal_edges.py', keepclusters[i], 1.5)

	trim <- system2('python3', args, stdout = TRUE, stderr = TRUE)
	trim <- grep('removing', trim, value=TRUE)
	if (length(trim) > 0) {
		for (j in 1:length(trim)) {
			if (!grepl('^removing', trim[j])) {
				trim[j] <- strsplit(trim[j], ';')[[1]][2]
			}
		}
		trim <- strsplit(trim, ' ')
		trim <- sapply(trim, function(x) x[2])
		geneReport$trimmed[geneReport$accession %in% trim] <- TRUE
	}
	
	
	########
	# BLAST: for those samples that fail some of the above tests, blast those samples against the others in this locus to see if blast confirms the problems. If it doesn't, then problem may be due to poor tree inference.
	
	geneFasta <- ungappedSeqFiles[grep(paste0('^', gene), basename(ungappedSeqFiles))]

	# create blast database
	dir.create('blastTemp')
	file.copy(geneFasta, paste0('blastTemp/', basename(geneFasta)))
	makeblastdb(file = paste0('blastTemp/', basename(geneFasta)), dbtype = 'nucl', args = '-parse_seqids')
	db <- blast(paste0('blastTemp/', basename(geneFasta)))

	# read fasta seq data
	seq <- readDNAStringSet(geneFasta)

	outlierTaxonInd <- which(geneReport$patristicOutlier == TRUE | geneReport$trimmed == TRUE)
	if (length(outlierTaxonInd)) {
		for (j in 1:length(outlierTaxonInd)) {
			# cat(j, ' ')
			
			taxRec <- tax[which(tax$accession == geneReport[outlierTaxonInd[j], 'accession']), ]
			
			acc <- taxRec[, 'accession']
			if (!acc %in% names(seq)) stop('focal accession not in fasta names.')
			
			# perform megablast
			blastsearch <- predict(db, seq[acc, ], BLAST_args = '-task dc-megablast')
			blastsearch$QueryID <- as.character(blastsearch$QueryID)
			blastsearch$SubjectID <- as.character(blastsearch$SubjectID)
			blastsearch <- blastsearch[order(blastsearch$Perc.Ident, decreasing = TRUE),]
			
			# remove blast to self
			blastsearch <- blastsearch[blastsearch$SubjectID != acc,]
			
			# add in taxonomy
			if (!all(blastsearch$SubjectID %in% tax$accession)) {
				#stop('not all blast accessions in tax table.')
				blastsearch <- blastsearch[blastsearch$SubjectID %in% tax$accession,]
			}
			blastsearch <- cbind.data.frame(blastsearch, do.call(rbind, lapply(blastsearch$SubjectID, function(x) tax[which(tax$accession == x), c('family','genus','species')])))
			
			# are the top 10 closest blast hits in the same family? in same genus?
			head(blastsearch)
			# what percent of the 10 closest hits are in the same family? 
			sameFamily <- tax[which(tax$family == taxRec[, 'family']), 'species']
			nFam <- 10
			if (length(sameFamily) < nFam) {
				# stop('nFam > number of taxa in family')	
				nFam <- length(sameFamily)
			}
			percentInFamily <- sum(blastsearch[1:nFam, 'species'] %in% sameFamily) / nFam
	
			sameGenus <- tax[which(tax$genus == taxRec[, 'genus']), 'species']
			nGenus <- 10
			if (length(sameGenus) < nGenus) {
				# stop('nGenus > number of taxa in genus')
				nGenus <- length(sameGenus)
			}
			percentInGenus <- sum(blastsearch[1:nGenus, 'species'] %in% sameGenus) / nGenus
			
			geneReport[outlierTaxonInd[j], 'blastFamilyPercent'] <- percentInFamily
			geneReport[outlierTaxonInd[j], 'blastGenusPercent'] <- percentInGenus
			
			rm(blastsearch)
		}
	}
		
	range(geneReport$blastFamilyPercent, na.rm = TRUE)
	range(geneReport$blastGenusPercent, na.rm = TRUE)
	table(geneReport$blastFamilyPercent < 0.8, useNA = 'always')
	table(geneReport$blastGenusPercent < 0.8, useNA = 'always')
	
	unlink('blastTemp', recursive = TRUE)
	
	
	
	########
	# ALIGNMENT-BASED METRICS
	# check that subspecies form exclusive clades
		
	aln <- readDNAStringSet(gsub('\\.tre|\\.treefile$', '.aln', keepclusters[i]))
	aln <- lapply(aln, as.character)	
	aln <- lapply(aln, function(x) strsplit(x, '')[[1]])
	
	# number of private nucleotides (nucleotides that don't show up elsewhere in that site in the alignment)
	seqMat <- do.call(rbind, aln)	
	nt <- c('A','T','G','C')
		
	percentCutoff <- 0.01
		
	nPrivateNT <- numeric(length(aln))
	names(nPrivateNT) <- names(aln)

	# for each nucleotide, find the sites that have that value in less than % individuals.
	for (j in 1:length(nt)) {
		seqMat2 <- seqMat == nt[j]
		sites <- which(colSums(seqMat2) / nrow(seqMat) < percentCutoff)
		ind <- rowSums(seqMat2[, sites])
		nPrivateNT <- nPrivateNT + ind
	}
	
	nPrivateNT <- nPrivateNT[order(names(nPrivateNT))]

	geneReport <- geneReport[order(geneReport$accession),]

	geneReport$percentPrivateAlleles <- nPrivateNT[geneReport$accession] / ncol(seqMat)
	
	# For duplicated taxa, identify the best accession to keep
	# we will not use percent private alleles because if one seq is much longer than the other, more private alleles in the long seq would be misleading
	# We will just focus on monophyly and 'trimmed' (overly long branch)
	dups <- tax$species[which(duplicated(tax$species) == TRUE)]
	if (length(dups) > 0) {
		for (j in 1:length(dups)) {
			ind <- which(geneReport$species == dups[j])
			tmp <- geneReport[ind, ]
			if (length(ind) > 2) stop()
			keep <- ind
			drop <- NA
			# in order of increasing priority
			if (any(tmp$patristicOutlier == TRUE & tmp$patristicOutlier == FALSE)) {
				drop <- which(tmp$patristicOutlier == TRUE)
			}
			if (any(tmp$familyMonophyletic == TRUE & tmp$familyMonophyletic == FALSE)) {
				drop <- which(tmp$familyMonophyletic == FALSE)
			}
			if (any(tmp$genusMonophyletic == TRUE & tmp$genusMonophyletic == FALSE)) {
				drop <- which(tmp$genusMonophyletic == TRUE)
			}
			if (any(tmp$trimmed == TRUE & tmp$trimmed == FALSE)) {
				drop <- which(tmp$trimmed == TRUE)
			}
			if (!is.na(drop)) stop()
			
			if (is.na(drop)) {
				# keep the longer seq
				len <- tax[which(tax$species == dups[j]), c('accession', 'seqLength')]
				if (length(unique(len$seqLength)) == 1) {
					drop <- len[2, 'accession']
				} else {
					drop <- len[which(len$seqLength != max(len$seqLength))[1], 'accession']
				}
			}
			
			geneReport[which(geneReport$accession == drop), 'dropAsDup'] <- TRUE
		}
		
		for (j in 1:length(dups)) {
			ind <- which(geneReport$species == dups[j])
			tmp <- geneReport[ind, ]
			if (!any(tmp$dropAsDup == TRUE)) stop()
			
		}
	}
	
	message('\t\t', length(which(geneReport$patristicOutlier == T)), ' taxa are patristic outliers.')
	message('\t\t', sum(sapply(split(geneReport$familyMonophyletic, geneReport$family), function(x) sum(x) == 0)), ' families are not monophyletic.')
	message('\t\t', sum(sapply(split(geneReport$genusMonophyletic, geneReport$genus), function(x) sum(x) == 0)), ' genera are not monophyletic.')
	message('\t\tRange in percent private alleles: ', round(range(geneReport$percentPrivateAlleles), 2)[1], ' - ', round(range(geneReport$percentPrivateAlleles), 2)[2])
	message('\t\tRange in terminal edge length: ', round(range(geneReport$terminalEdge), 6)[1], ' - ', round(range(geneReport$terminalEdge), 6)[2])
	
	reportList[[i]] <- geneReport

}





write.csv(do.call(rbind, reportList), 'Squamata_clusterTreeReport.csv', row.names=F)
	
	
# plot trees, highlighting outliers

# read in report
report <- read.csv('Squamata_clusterTreeReport.csv', stringsAsFactors=FALSE)

clusters <- sort(unique(report$locus))

for (i in 1:length(clusters)) {
	
	message('\t', clusters[i])
	
	tree <- read.tree(grep(paste0('^', clusters[i]), keepclusters, value=TRUE))
	
	# superCRUNCH may have added * to the end of some labels. Remove those. 
	tree$tip.label <- gsub('\\*', '', tree$tip.label)

	dat <- report[which(report$locus == clusters[i]),]

	clade <- 'Squamata'
	pdf(paste0('clusterTrees/', clusters[i], '.pdf'), width=8, height=Ntip(tree) * 0.05)

	# subset tonini tree to as many of these taxa as possible
	toniniSub <- keep.tip(tonini, intersect(tonini$tip.label, dat$species))
	
	# what is the basal most taxon?
	xx <- which(toniniSub$edge[,1] == Ntip(toniniSub) + 1)
	xx <- toniniSub$edge[xx, 2]

	if (length(geiger::tips(toniniSub, xx[1])) < length(geiger::tips(toniniSub, xx[2]))) {
		outgroup <- geiger::tips(toniniSub, xx[1])
	} else {
		outgroup <- geiger::tips(toniniSub, xx[2])
	}

	if (any(grepl('Sphenodon', dat$species))) {
		outgroup <- grep('Sphenodon', dat$species, value=T)
	}
	
	# reroot
	outgroup <- lapply(outgroup, function(x) dat[which(dat$species == x), 'accession'])
	outgroup <- unlist(outgroup)
	if (length(outgroup) == 1) {
		tree <- ladderize(root(tree, outgroup))
	} else {
		tree <- ladderize(root(tree, node=getMRCA(tree, outgroup)))
	}

	# reorder dat to match tree
	dat <- dat[sapply(tree$tip.label, function(x) which(dat$accession == x)),]

	# blue for family-level problems
	# red for genus-level problems
	# coral for subspecies-level problems

	edgeColor <- rep('black', nrow(dat))

	ind <- which(dat$trimmed == TRUE)
	edgeColor[ind] <- 'blue'

	ind <- which(dat$patristicOutlier == TRUE)
	tmp <- dat[ind, c('blastFamilyPercent', 'blastGenusPercent')]
	ind <- ind[which(tmp$blastFamilyPercent < 0.8 | tmp$blastGenusPercent < 0.5)] # unflag records that had good blast results
	edgeColor[ind] <- 'red'
	
	ind <- which(dat$trimmed == TRUE & dat$patristicOutlier == TRUE)
	edgeColor[ind] <- 'purple'
	
	table(edgeColor)
	
	tipFonts <- rep(1, Ntip(tree))
	tipFonts[which(edgeColor != 'black')] <- 2
	
	familyLabels <- paste0(dat$family, '/ ', dat$species, ' ', dat$accession)
	
	# for accessions that were marked for removal because they are duplicate taxa, label those too
	if (any(dat$dropAsDup == TRUE)) {
		familyLabels[which(dat$dropAsDup == TRUE)] <- paste0(familyLabels[which(dat$dropAsDup == TRUE)], '//dupDrop')
	}
	
	tree$tip.label <- familyLabels
	plot.phylo(tree, tip.color = edgeColor, font = tipFonts, cex = 0.25, edge.width=0.8, no.margin=TRUE, underscore=TRUE, label.offset=0.001)	
	
	dev.off()

}






# how many times does a particular sample come up?
head(report)

outliers <- report[which(report$patristicOutlier == TRUE),]

head(sort(table(outliers$species), decreasing=T))
inMoreThan10 <- sort(table(outliers$species), decreasing=T)
inMoreThan10 <- inMoreThan10[inMoreThan10 > 9]

for (i in 1:length(inMoreThan10)) {
	outliers[which(outliers$species == names(inMoreThan10)[i]),c('cluster','genbankID')]
}





