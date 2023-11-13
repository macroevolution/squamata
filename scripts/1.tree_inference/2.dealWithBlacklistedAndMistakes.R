# Given a set of accessions, determine whether we are currently using them, and if we are, what to do about it.

# 	These include:
# 		accessions blacklisted by Zaher et al. 2019
#		accessions that are known by Pyron to be mistakes in taxon identity

# This file contains all possible seqs from local genbank DB
taxonTableFile <- '~/Dropbox/squamatePhylo/2019/topologyScripts/pyphlawdMissingTaxa/Lepidosauria_8504.table'
taxonTable <- data.table::fread(taxonTableFile, data.table=FALSE)
colnames(taxonTable) <- c('V1', 'taxonID', 'locus', 'accession', 'taxon', 'definition', 'title')
head(taxonTable)

# Here, we will pull together all accessions in pyphlawd clusters that we are employing. 
clusterDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/fastaClusters_withNonNCBI'

# put together list of accessions descriptions per cluster
setwd(clusterDir)
clusters <- list.files(pattern = '\\.fa$')
clusterList <- vector('list', length(clusters))
names(clusterList) <- gsub('\\.fa', '', clusters)
for (i in 1:length(clusters)) {
	message('\t', clusters[i])
	seq <- ape::read.FASTA(clusters[i])
	clusterList[[i]] <- taxonTable[taxonTable$accession %in% names(seq), c('taxon', 'taxonID', 'accession', 'definition')]
}

head(lapply(clusterList, head))


# Read in accessions from supplementary table S6 of Zaher et al. 2019
zaherBlacklisted <- '~/Dropbox/squamatePhylo/2019/Zaher_2019_blacklistedAccessions.txt'
zaherBlacklisted <- scan(zaherBlacklisted, what = 'character', sep = '\n')

# fixable errors identified by Pyron
toniniFixable <- c('DQ990973', 'DQ990972', 'GQ895805', 'GQ895861', 'GQ895806', 'GQ895862', 'EU402840')
# which of these are we using?
intersect(toniniFixable, do.call(rbind, clusterList)$accession)

# DQ990972 is 16S, listed as Anelytropsis papillosus in our data, should be Xenosaurus platyceps
# EU402840 is RAG1, listed as Casarea dussumieri (a boa) in our data, it is associated with Causus_defilippi in Zaher et al. 
# BLAST returns similarity to a bunch of vipers. According to source publication, it should be Causus defilippi. 
# => Fix name in taxonTable, find replacement 
# => Causus defilippi does not currently have RAG1, so renaming is fine. 
# => Need a replacement RAG1 for Casarea. AY487408 should do it. 
## Repeating these action items below. 


toReplace <- matrix(nrow = 0, ncol = 4)
colnames(toReplace) <- c('accession', 'taxonID', 'taxon', 'cluster')

for (i in 1:length(zaherBlacklisted)) {
	
	acc <- zaherBlacklisted[i]
	clusterInd <- which(sapply(clusterList, function(x) acc %in% x$accession))
	if (length(clusterInd) > 0) {
		message(i)
		clusterName <- names(clusterList)[clusterInd]
		toReplace <- rbind(toReplace, c(clusterList[[clusterInd]][which(clusterList[[clusterInd]]$accession == acc), c('accession', 'taxonID', 'taxon')], clusterName))
		
	}
}

# This is a table of accessions to drop, and of accessions that could be used to replace
write.csv(toReplace, '~/Dropbox/squamatePhylo/2019/topologyScripts/ZaherBlacklisted_updates.csv', row.names=FALSE)

toReplace[toReplace[,1] %in% toniniFixable,]

# For these that were flagged, do any of them belong to genes that have multiple clusters, and is the species in the other cluster as well?

for (i in 1:nrow(toReplace)) {
	gene <- toReplace[i, 'cluster']
	gene <- gsub('(.+)_cluster\\d+$', '\\1', gene)
	if (length(grep(gene, names(clusterList))) > 1) {
		message('\tMultiple clusters for ', gene)
		otherClusters <- setdiff(grep(gene, names(clusterList), value=TRUE), toReplace[i, 'cluster'])
		if (toReplace[i, 'taxon'] %in% do.call(rbind, clusterList[otherClusters])$taxon == TRUE) {	
			message('\t\ttaxon found!')
		}	
	}
}

# none appear to be in several clusters


# We will now use the taxonTable to hopefully identify suitable replacements. If none are found, these accessions will need to be dropped.

for (i in 1:nrow(toReplace)) {
	
	zz <- taxonTable[which(taxonTable$taxon == toReplace[i, 'taxon']),]
	zz <- zz[which(!zz$accession %in% zaherBlacklisted),]
	zz
	gene <- gsub('(.+)_cluster\\d+$', '\\1', toReplace[i, 'cluster'])
	zz[union(grep(gene, zz$definition, ignore.case = TRUE), grep(gene, zz$accession, ignore.case = TRUE)),]
	toReplace[i,]
}

## for those taxa where there are a number of possible replacements, let's determine which seq is longest
library(reutils)
library(AnnotationBustR)

# Zamenis longissimus CYTB
acc <- c("HQ392548", "HQ392549", "HQ392550", "HQ392551", "HQ392552", "HQ392553", "HQ392554", "HQ392555", "HQ392556", "HQ392557", "HQ392558", "HQ392559", "HQ392560", "HQ392561", "HQ392562", "HQ392563", "HQ392564", "HQ392565", "HQ392566", "MH018691", "MH018693")
FindLongestSeq(acc)
# HQ392548 is longest.

# Zamenis longissimus ND4
acc <- c("KY495532", "KY495533", "KY495534", "KY495535", "KY495536", "KY495537", "KY495538", "KY495539", "KY495540", "KY495541", "KY495542", "KY495543", "KY495544", "KY495545", "KY495546", "KY495547", "KY495548", "KY495549", "KY495550", "KY495551", "KY495552", "KY495553", "KY495554", "KY495555", "KY495556", "KY495557", "KY495558", "KY495559")
FindLongestSeq(acc)
# KY495532 is longest

# Zamenis longissimus CMOS
acc <- c("KY495511", "KY495512", "KY495513", "KY495514", "KY495515", "KY495516", "KY495517", "KY495518", "KY495519", "KY495520", "KY495521", "KY495522", "KY495523", "KY495524", "KY495525")
FindLongestSeq(acc)


######################
# ACTION ITEMS
######################

# all blacklisted accessions will be dropped

# can replace Indotyphlops braminus RAG1 with GU902633
# can replace Dipsadoboa werneri CYTB with MH841943 (maybe new since GenBank download)
# can replace Manolepis putnami CYTB with JQ598936
# Conophis vittatus CMOS GQ895806 is a known error => Should be Conopsis nasus. Simply fix in taxonTable? <= ALREADY HANDLED
# can replace Clelia clelia CMOS with JQ598973
# can replace Manolepis putnami CYTB with JQ598936
# can replace Siphlophis compressus CYTB with SqCL_UMMZ_245073_CYTB
# Conopsis nasus CYTB GQ895861 is a known error => Should be Conophis vittatus. Simply fix in taxonTable? <= ALREADY HANDLED
# Conopsis nasus CMOS GQ895805 is a known error => Should be Conophis vittatus. Simply fix in taxonTable? <= ALREADY HANDLED
# Can replace Zamenis longissimus CYTB with HQ392548 
# Can replace Zamenis longissimus ND4 with KY495532
# Can replace Zamenis longissimus CMOS with KY495511

### Other action items based on Pyron's help:
# DQ990972 is 16S, listed as Anelytropsis papillosus in our data, should be Xenosaurus platyceps ALREADY HANDLED
# EU402840 is RAG1, listed as Casarea dussumieri (A boa) in our data, it is associated with Causus_defilippi in Zaher et al. 
# BLAST returns similarity to a bunch of vipers. According to source publication, it should be Causus defilippi. 
# => Fix name in taxonTable, find replacement 
# => Causus defilippi does not currently have RAG1, so renaming is fine. 
# => Need a replacement RAG1 for Casarea. AY487408 should do it. 



# PENDING: for Thamnodynastes pallidus SqCL_CHUNB52599, why do we only have ND4 and 16S? That taxon is blacklisted for CYTB and CMOS. If we had more SqCL loci, we could replace for CYTB. <= Due to poor data quality. 

# PENDING: Should we be using SqCL loci whenever possible, even if PyPhlawd selected a genbank accession over SqCL?
## No. SqCL mt genes are not necessarily great. 

# checking that those genes that I looked for replacements for, are not for taxa found in another cluster.
for (i in 1:nrow(toReplace)) {
	zz <- which(sapply(clusterList, function(x) toReplace[i, 'taxonID'] %in% x$taxonID) == T)
	genes <- gsub('(.+)_cluster\\d+$', '\\1', names(zz))
	message('\t', i, ': ', paste(genes, collapse = ' '))
}
# not a concern.



















