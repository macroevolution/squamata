# Search for third order minor clusters (Pyron's wording): Within a genus, are there any taxa or groups of taxa for which there is no overlap in loci?


alnDir <- '~/squam2020/fasta_March2020_alignments'
outfile <- '../locusOverlap_toDrop.csv'

setwd(alnDir)

require(Biostrings)

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


# first create a taxon x locus table
geneList <- vector('list', length(alnFiles))
names(geneList) <- gsub(clusterTag, '', alnFiles)

for (i in 1:length(alnFiles)) {	
	message('\t', alnFiles[i])
	geneList[[i]] <- names(readDNAStringSet(alnFiles[i]))
}

allSp <- sort(unique(taxonTable$repdbTaxon))

mat <- matrix(nrow = length(allSp), ncol = length(alnFiles))
colnames(mat) <- names(geneList)
rownames(mat) <- allSp

for (i in 1:length(geneList)) {
	mat[, i] <- allSp %in% taxonTable[which(taxonTable$locus == colnames(mat)[i]), 'repdbTaxon']
}

colSums(mat) # number of samples per locus
mean(rowSums(mat)) # average number of loci per taxon
# hist(rowSums(mat))

genera <- sapply(allSp, function(x) strsplit(x, '_|\\s+')[[1]][1])

# Now for each genus, determine which loci each taxon has, and identify groupings

res <- list()

for (i in 1:length(unique(genera))) {
	
	message('\t', i)
	if (length(which(genera == unique(genera)[i])) > 1) {
		
		tmp <- mat[which(genera == unique(genera)[i]),]
		
		# group taxa into non-overlapping locus coverage	
		locusGroups <- list()
		locusGroups[[1]] <- rownames(tmp)[1]
		remainingTaxa <- setdiff(rownames(tmp), locusGroups[[1]])
		
		for (j in 1:length(remainingTaxa)) {
			tmpTaxon <- remainingTaxa[j]
			overlaps <- FALSE
			for (k in 1:length(locusGroups)) {	
				if (is.vector(tmp[locusGroups[[k]],]) == 1) {
					colsums <- tmp[locusGroups[[k]],]
				} else {
					colsums <- colSums(tmp[locusGroups[[k]],])
				}
				if (length(intersect(which(colsums > 0), which(tmp[tmpTaxon,] == 1))) > 0) {
					locusGroups[[k]] <- c(locusGroups[[k]], tmpTaxon)
					overlaps <- TRUE
				}
			}
			if (!overlaps) {
				locusGroups[length(locusGroups) + 1] <- tmpTaxon
			}
		}
		
		# collapse groups that are linked together via a particular sample
		locusGroups <- lapply(locusGroups, unique)
		locusGroups <- lapply(locusGroups, sort)
		
		lapply(locusGroups, function(x) {
			if (is.vector(tmp[x,]) == 1) {
				colsums <- tmp[x,]
			} else {
				colsums <- colSums(tmp[x,])
			}
			sort(names(which(colsums > 0)))
		})
		
		#if (length(locusGroups) > 2) stop()
		
		for (j in 1:length(locusGroups)) {
			grp1 <- locusGroups[[j]]
			for (k in 1:length(locusGroups)) {
				if (j < k) {
					grp2 <- locusGroups[[k]]
					if (length(intersect(grp1, grp2)) > 0) {
						locusGroups[[j]] <- union(grp1, grp2)
						locusGroups[[k]] <- NA
					}
				}
			}
		}
		locusGroups <- locusGroups[!sapply(locusGroups, anyNA)]
		locusGroups <- locusGroups[order(lengths(locusGroups))]
		res[[i]] <- locusGroups
		names(res)[i] <- unique(genera)[i]
	}
}


# for each element in the list, there are groupings of taxa that overlap in loci. 
# if an element in the list has > 1 group, then there is some non-overlap in loci.
# In those cases, the groups are sorted, increasing by size, so if there are 2 non-overlapping groups, then the first item will be the smaller number of species that don't overlap with the majority. 
# zero groups means monotypic genus.

# How many groups per genus?
table(lengths(res))

# which genera?
qq <- which(lengths(res) > 1)
qq

# how many taxa don't overlap with the rest of their genus?
sum(sapply(qq, function(x) length(res[[x]][[1]])))

# which ones?
sapply(qq, function(x) res[[x]][[1]])

# example
lapply(res[['Hypsilurus']], function(x) {
	if (is.vector(mat[x,])) {
		colsums <- mat[x,]
	} else {
		colsums <- colSums(mat[x,])
	}
	sort(names(which(colsums > 0)))
})

# the smaller groups will be written to file to be dropped

toDrop <- sort(unlist(sapply(qq, function(x) res[[x]][[1]])))
names(toDrop) <- NULL
toDrop <- as.matrix(toDrop)

write.csv(toDrop, outfile, row.names = FALSE)

