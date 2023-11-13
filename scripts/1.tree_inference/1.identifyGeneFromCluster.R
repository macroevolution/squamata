# collect accession descriptions for a cluster and search for keywords within it to determine which gene it represents. 
# goal is, for each gene of interest, identify the cluster(s) that are best representatives for it.

require(ape)
require(pbapply)

clustersDir <- getwd()

clusters <- list.files(pattern='\\.fa$')

taxonTableFile <- list.files(gsub('/clusters', '', clustersDir), pattern = '\\.table', full.name=TRUE)
taxonTable <- data.table::fread(taxonTableFile, data.table=FALSE)
colnames(taxonTable) <- c('V1', 'taxonID', 'locus', 'accession', 'taxon', 'definition', 'title')

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

# read in bait gene names
baits <- list.files('~/Documents/pyphlawd/baited/baits/baitsLong', pattern='\\.fa$')
baits <- gsub('\\.fa', '', baits)

all(baits %in% names(geneList))

## switch here for SLC8A1/3

# put together list of accessions descriptions per cluster
clusterDesc <- vector('list', length(clusters))
names(clusterDesc) <- gsub('\\.fa', '', clusters)
for (i in 1:length(clusters)) {
	message('\t', clusters[i])
	seq <- read.FASTA(clusters[i])
	clusterDesc[[i]] <- taxonTable[taxonTable$accession %in% names(seq), 'definition']
}



# finding the percent of a cluster that has mention of the gene is helpful for identifying the clusters that are that baits in the majority 
baitClusters <- vector('list', length(baits))
names(baitClusters) <- baits
for (i in 1:length(baits)) {

	message('\t', baits[i])
	terms <- geneList[[baits[i]]]
	terms <- paste(terms, collapse='|')
	
	# for each cluster, identify percent of seqs that contain the bait name in the description
	xx <- sapply(clusterDesc, function(x) sum(grepl(terms, x, ignore.case=TRUE)) / length(x))
	
	# let's only consider >= 10% and also clusters that have at least 50 species
	baitsFound <- xx[which(xx > 0.1)]
	baitsFound <- baitsFound[which(lengths(clusterDesc[names(baitsFound)]) >= 50)]
	baitClusters[[i]] <- names(baitsFound)
	message('\t\tnumber of clusters found: ', length(baitsFound))
	message('\t\tpercentage of seqs matched: ', paste(round(baitsFound,2), collapse=' '))
	message('\t\tSizes: ', paste(lengths(clusterDesc[names(baitsFound)]), collapse=' '))

	# what do the seqs that do not contain the terms list?
	# we also leave out labels with "sequence added on" as those seem to not have any name.
	# clusterDesc[[names(baitsFound)[j]]][!grepl(paste0(terms, '|sequence added on'), clusterDesc[[names(baitsFound)[j]]], ignore.case=T)]

}

# Notes: 
# just a few baits match to a few clusters
# COI: 4 clusters
# CYTB: 2 clusters
# ND1: 4 clusters
# ND2: 4 clusters
# ND4: 3 clusters
# PRLR: 2 clusters

# move them for further investigation
destDir <- '~/Documents/pyphlawd/Squamata_baitClusters_percIdentity50/'
for (i in 1:length(baitClusters)) {
	filesToMove <- c(paste0(baitClusters[[i]], '.fa'), paste0(baitClusters[[i]], '.aln'))
	destFiles <- paste0(destDir, filesToMove)
	newNames <- paste0(destDir, names(baitClusters)[i], '_', filesToMove)
	file.copy(filesToMove, destFiles)
	file.rename(destFiles, newNames)
}


CFastaReader: Hyphens are invalid and will be ignored around line 2


	# ADNP
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.24
	# 	Sizes: 899
	# AHR
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.24
	# 	Sizes: 830
	# AKAP9
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.3
	# 	Sizes: 883
	# AMEL
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 1
	# 	Sizes: 138
	# BACH1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.31
	# 	Sizes: 815
	# BDNF
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.6
	# 	Sizes: 1777
	# BHLHB2
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.23
	# 	Sizes: 846
	# BMP2
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.35
	# 	Sizes: 1012
	# CAND1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.27
	# 	Sizes: 913
	# CARD4
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.27
	# 	Sizes: 800
	# CILP
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.25
	# 	Sizes: 854
	# cmos
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.79
	# 	Sizes: 3062
	# COI
	# 	number of clusters found: 4
	# 	percentage of seqs matched: 0.15 0.79 0.31 0.43
	# 	Sizes: 816 1982 1596 392
	# CXCR4
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.26
	# 	Sizes: 852
	# CYTB
	# 	number of clusters found: 2
	# 	percentage of seqs matched: 0.96 0.94
	# 	Sizes: 369 3071
	# DLL1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.25
	# 	Sizes: 869
	# ECEL
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.41
	# 	Sizes: 601
	# ENC1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.14
	# 	Sizes: 848
	# FSHR
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.3
	# 	Sizes: 952
	# FSTL5
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.24
	# 	Sizes: 876
	# GALR1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.23
	# 	Sizes: 877
	# GHSR
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.23
	# 	Sizes: 847
	# GPR37
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.25
	# 	Sizes: 858
	# HLCS
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.33
	# 	Sizes: 634
	# INHIBA
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.23
	# 	Sizes: 858
	# LRRN1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.92
	# 	Sizes: 262
	# LZTSS1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.21
	# 	Sizes: 822
	# MKL1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.29
	# 	Sizes: 868
	# MLL3
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.31
	# 	Sizes: 606
	# MSH6
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.26
	# 	Sizes: 884
	# ND1
	# 	number of clusters found: 4
	# 	percentage of seqs matched: 0.22 0.95 0.47 0.47
	# 	Sizes: 913 60 1596 392
	# ND2
	# 	number of clusters found: 4
	# 	percentage of seqs matched: 0.99 0.65 0.72 0.96
	# 	Sizes: 816 1596 404 392
	# ND4
	# 	number of clusters found: 3
	# 	percentage of seqs matched: 1 0.99 0.91
	# 	Sizes: 145 73 2437
	# NGFB
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.31
	# 	Sizes: 979
	# NKTR
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.71
	# 	Sizes: 507
	# NTF-3
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.66
	# 	Sizes: 1563
	# PDC
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.63
	# 	Sizes: 1151
	# PNN
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.41
	# 	Sizes: 962
	# PRLR
	# 	number of clusters found: 2
	# 	percentage of seqs matched: 0.97 1
	# 	Sizes: 918 61
	# PTGER4
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.36
	# 	Sizes: 907
	# PTPN
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.35
	# 	Sizes: 745
	# R35
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.53
	# 	Sizes: 1196
	# RAG1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.76
	# 	Sizes: 2969
	# RAG2
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.49
	# 	Sizes: 1235
	# rRNA_12S
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.91
	# 	Sizes: 3076
	# rRNA_16S
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.91
	# 	Sizes: 3627
	# SINCAIP
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.39
	# 	Sizes: 1028
	# SLC30A1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.27
	# 	Sizes: 925
	# SLC8A1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.11
	# 	Sizes: 923
	# SLC8A3
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.15
	# 	Sizes: 923
	# TRAF6
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.26
	# 	Sizes: 880
	# UBN1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.79
	# 	Sizes: 307
	# VCPIP1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.22
	# 	Sizes: 869
	# ZEB2
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.95
	# 	Sizes: 339
	# ZFP36L1
	# 	number of clusters found: 1
	# 	percentage of seqs matched: 0.26
	# 	Sizes: 873



## Specific to SLC8A1/3 paralogs -- ran pyphlawd with 70% identity threshold

setwd('~/Documents/pyphlawd/SLC8A/Lepidosauria_8504/clusters')
clustersDir <- getwd()

clusters <- list.files(pattern='\\.fa$')

taxonTableFile <- list.files(gsub('/clusters', '', clustersDir), pattern = '\\.table', full.name=TRUE)
taxonTable <- data.table::fread(taxonTableFile, data.table=FALSE)
colnames(taxonTable) <- c('V1', 'taxonID', 'locus', 'accession', 'taxon', 'definition', 'title')

# put together list of accessions descriptions per cluster
clusterDesc <- vector('list', length(clusters))
names(clusterDesc) <- gsub('\\.fa', '', clusters)
for (i in 1:length(clusters)) {
	message('\t', clusters[i])
	seq <- read.FASTA(clusters[i])
	clusterDesc[[i]] <- taxonTable[taxonTable$accession %in% names(seq), 'definition']
}



# finding the percent of a cluster that has mention of the gene is helpful for identifying the clusters that are that baits in the majority 
baits <- c('SLC8A1', 'SLC8A3')
baitClusters <- vector('list', length(baits))
names(baitClusters) <- baits
for (i in 1:length(baits)) {

	message('\t', baits[i])
	terms <- geneList[[baits[i]]]
	terms <- paste(terms, collapse='|')
	
	# for each cluster, identify percent of seqs that contain the bait name in the description
	xx <- sapply(clusterDesc, function(x) sum(grepl(terms, x, ignore.case=TRUE)) / length(x))
	
	# let's only consider >= 10% and also clusters that have at least 50 species
	baitsFound <- xx[which(xx > 0.1)]
	baitsFound <- baitsFound[which(lengths(clusterDesc[names(baitsFound)]) >= 50)]
	baitClusters[[i]] <- names(baitsFound)
	message('\t\tnumber of clusters found: ', length(baitsFound))
	message('\t\tpercentage of seqs matched: ', paste(round(baitsFound,2), collapse=' '))
	message('\t\tSizes: ', paste(lengths(clusterDesc[names(baitsFound)]), collapse=' '))

	# what do the seqs that do not contain the terms list?
	# we also leave out labels with "sequence added on" as those seem to not have any name.
	# clusterDesc[[names(baitsFound)[j]]][!grepl(paste0(terms, '|sequence added on'), clusterDesc[[names(baitsFound)[j]]], ignore.case=T)]

}

	SLC8A1
		number of clusters found: 1
		percentage of seqs matched: 0.27
		Sizes: 923
	SLC8A3
		number of clusters found: 1
		percentage of seqs matched: 0.25
		Sizes: 896

> baitClusters
$SLC8A1
[1] "cluster13871"

$SLC8A3
[1] "cluster14555"



