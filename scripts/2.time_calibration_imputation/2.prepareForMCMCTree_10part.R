## Generate input files for MCMCTREE

# MCMCTREE
#- write new tree file encoded with fossil calibrations

# devtools::install_github('PuttickMacroevolution/MCMCtreeR')

require(ape)
require(MCMCtreeR)

setwd('~/Dropbox/squamatePhylo/2019')

## -----------------------------------------------------
# INPUTS
seqFile <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/SortaDate/sorta_concatenated_v5/concat_n10.phy'
partitionsFile <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/SortaDate/sorta_concatenated_v5/concat_n10.partitions'
treefile <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/tier1/fulltree/bestTrees/fulltree_default_con_1_raxmlOpt.raxml.bestTree'
genomDatFile <- 'squamate_phylogenomics_v11.csv'
taxonTableFile <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/sqCL_taxa_repDBmatching.csv'
higherTaxonomyFile <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/fasta_March2020/metatable3.csv'
fossilfile <- 'squamate_fossil_calibrations_v5_20230522_MEHJ.csv' # response to peer review

## -----------------------------------------------------
# OUTPUTS
mcmctreeDir <- '~/Dropbox/squamatePhylo/2019/part10/withOutgroups/ddbd_ARS3_June2023'

## ------------------------------------------------------

taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
# manual matching supercedes repDB auto matching, unless no manual match found
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'repDBtaxon']
head(taxonTable)
table(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch))
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'species']

higherTaxonomy <- read.csv(higherTaxonomyFile)
higherTaxonomy[which(higherTaxonomy$genus %in% c('Myersophis', 'Oxyrhabdium', 'Hologerrhum', 'Cyclocorus' )), 'family'] <- 'Lamprophiidae'
higherTaxonomy[which(higherTaxonomy$genus == 'Xylophis'), 'family'] <- 'Pareidae'



# Read in seq and tree files, and prune/rename to common set of genus_species
genom <- read.tree(treefile)
#dat <- read.csv(genomDatFile, stringsAsFactors=FALSE)

# partitions file
part <- scan(partitionsFile, what = 'character', sep = '\n')
part <- gsub('(\\d)(-)(\\d)', '\\1|\\3', part)
part <- lapply(part, function(x) strsplit(x, split = '=|\\|')[[1]])
names(part) <- sapply(part, function(x) x[1])
part <- lapply(part, function(x) x[2:length(x)])
part <- lapply(part, as.numeric)
part

# sequence data: process and get genus_species from sample names
seq <- scan(seqFile, what = 'character', sep = '\n')
# first line is number of species/sites
spBySites <- seq[[1]]
spBySites <- strsplit(spBySites, split = '\t')[[1]]
seq <- seq[2:length(seq)]

# split out sample names
seqList <- vector('list', length = length(seq))

for (i in 1:length(seq)) {

	tmp <- strsplit(seq[[i]], split = '\t')[[1]]
	seqList[[i]] <- tmp[[2]]
	names(seqList)[i] <- tmp[[1]]
}

# remove outgroups
genomOutgroups <- c('taeGut2', 'ISIS373002', 'UMFS-10956c', 'allMis1', 'chrPic1', 'galGal5', 'hg38', 'H20145a')
genomOutgroups %in% names(seqList)
genomOutgroups <- intersect(genomOutgroups, names(seqList))

# ##########################

sapply(genomOutgroups, function(x) taxonTable[which(taxonTable$sample == x)[1], 'manualMatch'])

# keep Homo sapiens, Gallus gallus and Alligator as outgroups: hg38 and galGal5, allMis1

seqList <- seqList[setdiff(names(seqList), setdiff(genomOutgroups, c('hg38', 'galGal5', 'allMis1')))]

table(names(seqList) %in% taxonTable$sample)

# drop the few that are not matching
seqList <- seqList[intersect(names(seqList), taxonTable$sample)]

newLabels <- sapply(names(seqList), function(x) taxonTable[which(taxonTable$sample == x)[1], 'manualMatch'])
newLabels <- gsub('\\s+', '_', newLabels)
table(newLabels == '' | is.na(newLabels))

names(seqList) <- newLabels

# drop Rhineura
seqList <- seqList[setdiff(names(seqList),'Rhineura_floridana')]


# prune phylogeny to the same set of taxa

# how many seq species are not in the tree?
length(setdiff(names(seqList), genom$tip.label))
length(setdiff(genom$tip.label, names(seqList)))
length(intersect(genom$tip.label, names(seqList)))

commonTaxa <- intersect(genom$tip.label, names(seqList))

# keep only common taxa in seq data
seqList <- seqList[c(commonTaxa, 'Homo_sapiens', 'Gallus_gallus', 'Alligator_mississippiensis')]

# Root tree and prune
genom <- root(genom, 'Sphenodon_punctatus', resolve.root=TRUE)
# plot(ladderize(genom), no.margin = TRUE, cex=0.5)
subtree <- keep.tip(genom, commonTaxa)

#plot(ladderize(subtree), no.margin=TRUE, cex=0.3)

# check that both seq and tree have identical taxa
identical(sort(names(seqList)), sort(subtree$tip.label))
setdiff(names(seqList), subtree$tip.label)
setdiff(subtree$tip.label, names(seqList))

# graft on outgroups
outgroup <- read.tree(text = '(((Gallus_gallus, Alligator_mississippiensis), squamata), Homo_sapiens);')
#outgroup <- read.tree(text = '((Gallus_gallus, Alligator_mississippiensis), Homo_sapiens);')
#subtree <- bind.tree(x = subtree, y = outgroup, where = 'root', position = 0.2)
subtree <- bind.tree(x = outgroup, y = subtree, where = which(outgroup$tip.label == 'squamata'))
#subtree <- root(subtree, outgroup = 'Homo_sapiens', resolve.root = TRUE)
#subtree <- root(subtree, outgroup = 'outgroup')
#subtree <- drop.tip(subtree, 'outgroup')
#subtree <- root(subtree, outgroup = 'Homo_sapiens')
is.binary(subtree)

subtree <- read.tree(text = write.tree(ladderize(subtree)))

plot(ladderize(subtree), no.margin=TRUE, cex=0.3)

##------------------------------------------------
# write seqs for separate sanity check tree inference

# # treeInfDir <- '~/Dropbox/squamatePhylo/2019/part10/withOutgroups/treeInferenceTest/'

# for (i in 1:length(seqList)) {
	# write(paste0('>', names(seqList)[i]), file = paste0(treeInfDir, 'seqs.aln'), append = ifelse(i == 1, FALSE, TRUE))
	# write(seqList[[i]], file = paste0(treeInfDir, 'seqs.aln'), append = TRUE)
# }

# # partition file
# write(paste0('DNA, locus1 = ', part[[1]][1], '-', part[[1]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = FALSE)
# write(paste0('DNA, locus2 = ', part[[2]][1], '-', part[[2]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = TRUE)
# write(paste0('DNA, locus3 = ', part[[3]][1], '-', part[[3]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = TRUE)
# write(paste0('DNA, locus4 = ', part[[4]][1], '-', part[[4]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = TRUE)
# write(paste0('DNA, locus5 = ', part[[5]][1], '-', part[[5]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = TRUE)
# write(paste0('DNA, locus6 = ', part[[6]][1], '-', part[[6]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = TRUE)
# write(paste0('DNA, locus7 = ', part[[7]][1], '-', part[[7]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = TRUE)
# write(paste0('DNA, locus8 = ', part[[8]][1], '-', part[[8]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = TRUE)
# write(paste0('DNA, locus9 = ', part[[9]][1], '-', part[[9]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = TRUE)
# write(paste0('DNA, locus10 = ', part[[10]][1], '-', part[[10]][2]), file = paste0(treeInfDir, 'partitions.iqtree'), append = TRUE)

# # write topology
# write.tree(unroot(subtree), file = paste0(treeInfDir, 'topology.tre'))

# ~/iqtree2 -s seqs.aln -spp partitions.iqtree -m GTR+G -nt 2 -te topology.tre --runs 3 -pre part10 -o Homo_sapiens

## -----------------------------------------------
# Run CorrTest (assume tree was inferred in previous step)
setwd(treeInfDir)

source('~/CorrTest/code/rate.CorrTest.R')
infTree <- read.tree('part10.treefile')
outgroupTips <- gsub('\\s+', '_', sapply(genomOutgroups, function(x) taxonTable[which(taxonTable$sample == x)[1], 'manualMatch']))
outgroupTips <- intersect(outgroupTips, infTree$tip.label)
infTree <- root(infTree, 'Homo_sapiens')
# can't have polytomies, drop all but one outgroup
#infTree <- drop.tip(infTree, outgroupTips[2:length(outgroupTips)])
plot(ladderize(infTree), cex = 0.5)
rate.CorrTest(brlen_tree = infTree, outgroup = intersect(infTree$tip.label, outgroupTips), sister.resample = 0, outputFile = "CorrTest.txt")

# A significant p-value indicates that we reject the null hypothesis that rates are independent. 

## -----------------------------------------------
setwd('~/Dropbox/squamatePhylo/2019')



# getSpanningTips
#	returns a pair of tips that span a given node. 
#	if the node is terminal, includes "NA"
getSpanningTips <- function(phy, node){	
	if (node <= length(phy$tip.label)){
		return(c(phy$tip.label[node], 'NA'));
	}else{
		dnodes <- phy$edge[,2][phy$edge[,1] == node];

		while (dnodes[1] > length(phy$tip.label)){
			dnodes[1] <- phy$edge[,2][phy$edge[,1] == dnodes[1]][1];
		}
		while (dnodes[2] > length(phy$tip.label)){
			dnodes[2] <- phy$edge[,2][phy$edge[,1] == dnodes[2]][1];
		}		
		dset <- phy$tip.label[dnodes];
		return(dset);
	}	
}

# new tree file with no branch lengths, and fixed root age
subtree$node.label <- NULL
subtree2 <- subtree

# remove branch lengths
subtree$edge.length <- NULL



## -----------------------------------------------
# Generate MCMCTREE output files


# For seq file, MCMCTREE requires that different loci follow each other in file, rather than be concatenated.

# write new phylip file for use with PAML mcmctree
newSeqFile <- paste0(mcmctreeDir, '/part10_renamed.phy')

for (i in 1:length(part)) {
	
	# create new header lines
	# first line is # of species, # of sites, and option G to denote parititoned analysis
	header1 <- paste(length(seqList), (sapply(part, diff) + 1)[i], sep = ' ')

	# write new file
	write(header1, file = newSeqFile, append = ifelse(i == 1, FALSE, TRUE))

	for (j in 1:length(seqList)) {
		write(names(seqList)[j], file = newSeqFile, append = TRUE)
		write(substring(seqList[[j]], part[[i]][1], part[[i]][2]), file = newSeqFile, append = TRUE)
	}
	
	# add a blank line between loci
	write('', file = newSeqFile, append = TRUE)
}



fossildat <- read.csv(fossilfile, stringsAsFactors=FALSE)

fossildat <- fossildat[which(fossildat[,1] != ''),]
fossildat <- fossildat[which(fossildat$use != ''),]
fossildat <- fossildat[which(fossildat$use == 'yes'),]

# apply changes for revision
toDrop <- which(fossildat$revision == 'no')
fossildat <- fossildat[- toDrop,]

# subset table to only fossils for which spanning taxa are in tree
useFossils <- c()

for (i in 1:nrow(fossildat)) {
	spanningTaxa <- strsplit(fossildat[i, 'spanning.taxa'], split = '\\|')[[1]]
	if (all(spanningTaxa %in% subtree$tip.label)) {
		useFossils <- c(useFossils, i)
	}
}

# not in tree, but otherwise flagged for use
fossilsNotUsed <- fossildat[setdiff(1:nrow(fossildat), useFossils), ]

fossildat <- fossildat[useFossils, ]

# identify constrained nodes for plotting purposes
fossilNodes <- c()
for (i in 1:nrow(fossildat)) {
	spanningTaxa <- strsplit(fossildat[i, 'spanning.taxa'], split = '\\|')[[1]]
	fossilNodes <- c(fossilNodes, getMRCA(subtree, spanningTaxa))
}

plot(ladderize(subtree2), no.margin=T, cex=0.5)
nodelabels(text = rep('', length(fossilNodes)), node = fossilNodes, frame='circle', cex=0.3, bg = 'orange')

subtree3 <- subtree2
famLabels <- sapply(subtree3$tip.label, function(x) higherTaxonomy[which(higherTaxonomy$genus == strsplit(x, '_')[[1]][1])[1], 'family'])
subtree3$tip.label <- paste0(famLabels, '/ ', subtree3$tip.label)
# pdf('fossilPlacements.pdf', width=8, height=15)
# plot(ladderize(subtree3), no.margin=T, cex=0.5)
# nodelabels(text = rep('', length(fossilNodes)), node = fossilNodes, frame='circle', cex=0.3, bg = 'orange')
# nodelabels(text = fossildat$nodeLabel, node = fossilNodes, adj = c(1.3,1.2), frame='none', cex=0.5, bg = 'transparent', font=2, col='blue')
# dev.off()


# remove branch lengths from phylo object
subtree$edge.length <- NULL

# subtree <- root(subtree, 'Homo_sapiens', resolve.root = TRUE)
is.binary(subtree)


# When we have both min and max fossil bounds, we will apply uniform priors with soft bounds for min/max fossils (lower bound 1%, upper bound 5%)

# When we lack an upper bound, we will skew normal
# with soft minimum bound of 0.1% and soft upper bound of 95%




# the MCMCTreeR function seems to be optimized to work with all calibrations in one go, producing a properly calibrated tree at the end. We will therefore prepare different vectors for the different parameters. 


# ----------------------------
## STRATEGY S3 -- skew-normal distribution

fossilOutput <- paste0(mcmctreeDir, '/fossilCalibrations_S3.tre')

fossilParams <- list()

for (i in 1:nrow(fossildat)) {
	
	cat(i, '\t')
	
	spanningTaxa <- strsplit(fossildat[i, 'spanning.taxa'], split = '\\|')[[1]]
			
	minAge <- as.numeric(fossildat[i, 'min.age'])/100 # changing to units of 100 mya
	maxAge <- as.numeric(fossildat[i, 'max.age'])/100
	
	outputFiles <- FALSE
	# if adding the last calibration, have the function write to file
	if (i == nrow(fossildat)) {
		outputFiles <- TRUE
	} else {
		outputFiles <- FALSE
	}
	
	if (i == 1) {
		# create initial object, that will then be added to
		if (is.na(maxAge)) {

			fossilRes <- estimateSkewNormal(minAge = minAge, maxAge = minAge+0.2, phy = subtree, monoGroups = spanningTaxa, shape = 7, scale = 0.05, addMode = 0, maxProb = 0.95, minProb = 0.001, estimateScale = TRUE, estimateShape = FALSE, estimateMode = FALSE, plot=F, writeMCMCtree = FALSE)
		
		} else {
			fossilRes <- estimateBound(minAge = minAge, maxAge = maxAge, minProb = 0.01, phy = subtree, rightTail = 0.05, monoGroups = spanningTaxa, writeMCMCtree = FALSE, plot = FALSE)
		}
		
	} else {
		
		# add to the initial data object
		if (is.na(maxAge)) {
			
			fossilRes <- estimateSkewNormal(minAge = minAge, maxAge = minAge+0.2, phy = fossilRes$apePhy, monoGroups = spanningTaxa, shape = 7, scale = 0.05, addMode = 0, maxProb = 0.95, minProb = 0.001, estimateScale = TRUE, estimateShape = FALSE, estimateMode = FALSE, plot=F, writeMCMCtree = outputFiles, MCMCtreeName = fossilOutput)
			 			 
		} else {
			fossilRes <- estimateBound(minAge = minAge, maxAge = maxAge, minProb = 0.01, phy = fossilRes$apePhy, rightTail = 0.05, monoGroups = spanningTaxa, writeMCMCtree = outputFiles, plot = F, MCMCtreeName = fossilOutput)
		}			
	}
	fossilParams[[i]] <- fossilRes$parameters
}


