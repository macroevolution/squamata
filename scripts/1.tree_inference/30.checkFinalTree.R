# # Examine finalized tree
# # 	Check for monophyly of families and genera
# # 	Check for monophyly of subspecies with their nominal taxon


require(ape)
require(phytools)
require(phangorn)
require(MonoPhy)
require(pbapply)

if (grepl('pascal', getwd())) {
	basedir <- '~/Dropbox/Oz_Crown_Ages/'
} else {
	basedir <- '~/Dropbox/squamate_tree/'
}


familiesOutDir <- paste0(basedir, 'tier1/families/')
alnDir <- paste0(basedir, 'finalAlignments_June2020/renamed')
alnFiles <- list.files(alnDir, pattern = '\\.aln$', full.names = FALSE)

fullAln <- paste0(basedir, 'finalAlignments_June2020/renamed/concatenated.phy')

setwd(basedir)

# Read in genomic constraint tree and rename tips from samples to taxon names
genomTree <- read.tree(paste0(basedir, 'data/all_bootstraps/all.concat_ind0.01_loci0.05_all_n5185.constraint.tre'))
genomTree$node.label <- NULL

# remove outgroups
genomOutgroups <- c('taeGut2', 'ISIS373002', 'UMFS-10956c', 'allMis1', 'chrPic1', 'galGal5', 'hg38', 'H20145a')
genomOutgroups %in% genomTree$tip.label
genomOutgroups <- genomOutgroups[genomOutgroups %in% genomTree$tip.label]
genomTree <- root(genomTree, outgroup = genomOutgroups)
genomTree <- ladderize(genomTree)
genomTree <- drop.tip(genomTree, genomOutgroups)


taxonTable <- read.csv('~/Dropbox/Oz_Crown_Ages/sqCL_taxa_repDBmatching.csv', stringsAsFactors=FALSE)
# manual matching supercedes repDB auto matching, unless no manual match found
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'repDBtaxon']
head(taxonTable)
table(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch))
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'species']

table(genomTree$tip.label %in% taxonTable$sample)
genomTree$tip.label[!genomTree$tip.label %in% taxonTable$sample]

newLabels <- sapply(genomTree$tip.label, function(x) taxonTable[which(taxonTable$sample == x)[1], 'manualMatch'])
table(newLabels == '' | is.na(newLabels))
newLabels <- gsub('\\s+', '_', newLabels)
names(newLabels) <- NULL
genomTree$tip.label <- newLabels

# Drop a few taxa from genomic constraint
genomTree <- drop.tip(genomTree, c('Bothrocophias_hyoprora', 'Corytophanes_hernandesii', 'Lophognathus_gilberti', 'Rhineura_floridana'))



# read in table that allows us to interpret tip labels
taxonTableFile <- 'fasta_March2020/metatable2.csv'
if (file.exists(taxonTableFile)) {
	taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
} else {
	taxonTableFile <- '~/squam2020/metatable2.csv'
	taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
}
taxonTable$repdbTaxon <- gsub('\\s+', '_', taxonTable$repdbTaxon)

# Corrections
taxonTable[which(taxonTable$genus %in% c('Myersophis', 'Oxyrhabdium', 'Hologerrhum', 'Cyclocorus' )), 'family'] <- 'Lamprophiidae'
taxonTable[which(taxonTable$genus == 'Xylophis'), 'family'] <- 'Pareidae'

excludeTaxa <- sort(unique(c(taxonTable[which(taxonTable$genus %in% c('Micrelaps', 'Buhoma')), 'repdbTaxon'], taxonTable[which(taxonTable$family %in% c('Dibamidae')), 'repdbTaxon'])))

families <- sort(unique(taxonTable$family))

# -----------------------------------------
## READ IN TREES WITH BEST LIKELIHOODS

setwd(basedir)

contree <- read.tree('tier1/fulltree/constrained/fulltreeInferences/fulltree_default_con_1_raxmlOpt.raxml.bestTree')
contree <- root(contree, 'Sphenodon_punctatus', resolve.root = TRUE)
contree <- ladderize(contree)

uncontree <- read.tree(paste0(basedir, 'tier1/fulltree/unconstrained/fulltreeInferences/fulltree_default_uncon_raxml_1_raxmlOpt.raxml.bestTree'))
uncontree <- root(uncontree, 'Sphenodon_punctatus', resolve.root = TRUE)
uncontree <- ladderize(uncontree)


## Did the family topologies remain fixed as expected?

bestFamilyDir <- paste0(basedir, 'tier1/fulltree/bestFamilyTrees')

## CONSTRAINED

# read in best-estimate clade trees
clades <- list.files(bestFamilyDir, pattern = '\\.tre')
clades <- sort(unique(gsub('_con|_uncon|\\.tre', '', clades)))

cladeList <- vector('list', length(clades))
names(cladeList) <- clades

famFiles <- list.files(bestFamilyDir, full.name = TRUE)
for (i in 1:length(clades)) {
	tmp <- famFiles[grep(clades[i], basename(famFiles))]
	if (length(tmp) == 2) {
		tmp <- tmp[grep('_con', basename(tmp))]
	} else if (length(tmp) > 2) stop()
	
	tmp <- read.tree(tmp)	

	tax <- unique(taxonTable[which(taxonTable$family == clades[i]), 'repdbTaxon'])
	
	toDrop <- setdiff(tmp$tip.label, tax)
	tmp <- drop.tip(tmp, toDrop)

	cladeList[[i]] <- tmp
}

for (i in 1:length(cladeList)) {
	
	#check <- all.equal.phylo(extract.clade(contree, getMRCA(contree, cladeList[[i]]$tip.label)), cladeList[[i]], use.edge.length = FALSE)
	check <- all.equal.phylo(keep.tip(contree, cladeList[[i]]$tip.label), cladeList[[i]], use.edge.length = FALSE)
	if (!check) {
		message('\t' ,i, ' ', clades[i], ' mismatch!')
		# famCo <- cophylo(keep.tip(contree, cladeList[[i]]$tip.label), cladeList[[i]])
		# plot(famCo)
	}
}

## Opluridae and Trogonophidae did not pass check, but that's because there was no tip in the genomic constraint to hang them on, so they were left out. Small number of tips.

## UNCONSTRAINED

# read in best-estimate clade trees
clades <- list.files(bestFamilyDir, pattern = '\\.tre')
clades <- sort(unique(gsub('_con|_uncon|\\.tre', '', clades)))

cladeList <- vector('list', length(clades))
names(cladeList) <- clades

famFiles <- list.files(bestFamilyDir, full.name = TRUE)
for (i in 1:length(clades)) {
	tmp <- famFiles[grep(clades[i], basename(famFiles))]
	if (length(tmp) == 2) {
		tmp <- tmp[grep('_uncon', basename(tmp))]
	} else if (length(tmp) > 2) stop()
	
	tmp <- read.tree(tmp)	

	tax <- unique(taxonTable[which(taxonTable$family == clades[i]), 'repdbTaxon'])
	
	toDrop <- setdiff(tmp$tip.label, tax)
	tmp <- drop.tip(tmp, toDrop)

	cladeList[[i]] <- tmp
}

for (i in 1:length(cladeList)) {
	
	#check <- all.equal.phylo(extract.clade(uncontree, getMRCA(contree, cladeList[[i]]$tip.label)), cladeList[[i]], use.edge.length = FALSE)
	check <- all.equal.phylo(keep.tip(uncontree, cladeList[[i]]$tip.label), cladeList[[i]], use.edge.length = FALSE)
	if (!check) {
		message('\t' ,i, ' ', clades[i], ' mismatch!')
		# famCo <- cophylo(keep.tip(contree, cladeList[[i]]$tip.label), cladeList[[i]])
		# plot(famCo)
	}
}

# all passed!


# -------------------------------------------------------

# What is the difference between the best constrained tree, in terms of IQTREE and re-optimized raxml-ng?

contree_iq <- read.tree('tier1/fulltree/constrained/fulltreeInferences/fulltree_default_con_1.treefile')
contree_iq <- root(contree_iq, 'Sphenodon_punctatus', resolve.root = TRUE)
contree_iq <- ladderize(contree_iq)

all.equal.phylo(contree, contree_iq, use.edge.length = FALSE) # TRUE

identical(contree$tip.label, contree_iq$tip.label) # TRUE

plot(contree_iq$edge.length, contree$edge.length, xlab = 'iqtree edge lengths', ylab = 'raxml reoptimized', cex =0.75, asp = 1, main = 'branch length comparison')
abline(a = 0, b = 1, lty = 2)

sum(contree_iq$edge.length)
sum(contree$edge.length)

# raxml-ng re-optimized has slightly less branch length.




# ---------------------------------------------------
# Exploration of monophyly

tax <- cbind(	species = contree$tip.label,
				genus = sapply(contree$tip.label, function(x) taxonTable[which(taxonTable$repdbTaxon == x)[1], 'genus']),
				family = sapply(contree$tip.label, function(x) taxonTable[which(taxonTable$repdbTaxon == x)[1], 'family']))

monoAssess <- AssessMonophyly(contree, tax[contree$tip.label,])

# genus monophyly

pdf('tier1/fulltree/constrained/bestConstrainedTree_genusMonophy.pdf', width = 15, height = 260)

PlotMonophyly(monoAssess, contree, taxlevels = 'genus', plot.type = 'monophyly', monocoll = TRUE, cex = 0.8, label.offset = 0.01, no.margin = TRUE, adj.names=0)
title(main = 'constrained tree', line = -10, cex.main = 2)

legend('topleft', legend = c('monophyletic', 'non-monophyletic', 'intruder/outlier'), fill = c('#5aae61', '#c2a5cf', '#762a83'))
dev.off()

# family monophyly

pdf('tier1/fulltree/constrained/bestConstrainedTree_familyMonophy.pdf', width = 15, height = 65)

PlotMonophyly(monoAssess, contree, taxlevels = 'family', plot.type = 'monophyly', monocoll = TRUE, cex = 0.8, label.offset = 0.01, no.margin = TRUE, adj.names=0)
title(main = 'constrained tree', line = -10, cex.main = 2)

legend('topleft', legend = c('monophyletic', 'non-monophyletic', 'intruder/outlier'), fill = c('#5aae61', '#c2a5cf', '#762a83'))
dev.off()


# ---------------------------------------------------
# Plot full tree, coloring according to family and genus level monophyly

# Also plot branches that represent genomic constraint with a thicker line

families <- unique(tax[, 'family'])
families <- setdiff(families, 'Sphenodontidae')

genera <- unique(tax[,'genus'])
genera <- setdiff(genera, 'Sphenodon')

contree2 <- contree
contree2 <- drop.tip(contree2, 'Sphenodon_punctatus')

familyGenus <- data.frame(tip = contree2$tip.label, family = tax[contree2$tip.label, 'family'], genus = tax[contree2$tip.label, 'genus'], familyMono = NA, genusMono = NA)

for (i in 1:length(families)) {
	spInGroup <- tax[which(tax[, 'family'] == families[i]), 'species']
	if (length(spInGroup) > 1) {
		tmp <- extract.clade(contree2, getMRCA(contree2, spInGroup))
		if (length(spInGroup) == Ntip(tmp)) {
			familyGenus[which(familyGenus$tip %in% spInGroup), 'familyMono'] <- 'TRUE'
		} else {
			message('\t', i, ' ', families[i], ' not monophyletic.')
			familyGenus[which(familyGenus$tip %in% spInGroup), 'familyMono'] <- 'FALSE'
		}
	} else {
		familyGenus[which(familyGenus$tip %in% spInGroup), 'familyMono'] <- 'TRUE'
	}
}

# genera
for (i in 1:length(genera)) {
	spInGroup <- tax[which(tax[, 'genus'] == genera[i]), 'species']
	if (length(spInGroup) > 1) {
		tmp <- extract.clade(contree2, getMRCA(contree2, spInGroup))
		if (length(spInGroup) == Ntip(tmp)) {
			familyGenus[which(familyGenus$tip %in% spInGroup), 'genusMono'] <- 'TRUE'
		} else {
			message('\t', i, ' ', genera[i], ' not monophyletic.')
			familyGenus[which(familyGenus$tip %in% spInGroup), 'genusMono'] <- 'FALSE'
		}
	} else {
		familyGenus[which(familyGenus$tip %in% spInGroup), 'genusMono'] <- 'TRUE'
	}
}

table(familyGenus$familyMono)
table(familyGenus$genusMono)

familyGenusCon <- familyGenus

tipColor <- rep('black', Ntip(contree2))

# lack of genus monophyly -> orange
# lack of family monophyly -> blue
# lack of both genus and family monophyly -> purple

tipColor[which(familyGenusCon$genusMono == FALSE & familyGenusCon$familyMono == TRUE)] <- 'orange'
tipColor[which(familyGenusCon$genusMono == TRUE & familyGenusCon$familyMono == FALSE)] <- 'blue'
tipColor[which(familyGenusCon$genusMono == FALSE & familyGenusCon$familyMono == FALSE)] <- 'purple'

table(tipColor)

# identify branches that are consistent with constraint tree
genomOverlap <- intersect(contree2$tip.label, genomTree$tip.label)

edgeVec <- pbsapply(contree2$edge[,2], function(x) any(genomOverlap %in% geiger::tips(contree2, x)), cl = 4)
edgeVec <- as.character(edgeVec)
table(edgeVec)
edgeVec[which(edgeVec == 'TRUE')] <- 'darkorchid3'
edgeVec[which(edgeVec == 'FALSE')] <- 'black'
table(edgeVec)

edgeLwd <- edgeVec
edgeLwd[which(edgeLwd == 'darkorchid3')] <- 3
edgeLwd[which(edgeLwd == 'black')] <- 1
edgeLwd <- as.numeric(edgeLwd)
table(edgeLwd)

newLabels <- contree2$tip.label
for (i in 1:length(newLabels)) {
		
	newLabels[i] <- paste0(tax[contree2$tip.label[i], 'family'], '_>_', newLabels[i])
	if (contree2$tip.label[i] %in% genomTree$tip.label) {
		newLabels[i] <- paste0(newLabels[i], ' ***')
	}
}

contree2$tip.label <- newLabels

pdf('tier1/fulltree/bestConstrainedTree_monophylyColors.pdf', width = 20, height = 700)

plot.phylo(contree2, cex = 0.5, no.margin = TRUE, label.offset = 0.001, tip.color = tipColor, edge.color = edgeVec, edge.width = edgeLwd)

dev.off()

pdf('tier1/fulltree/bestConstrainedTree_showConstraints.pdf', width = 100, height = 100)

plot.phylo(contree2, type = 'fan', show.tip.label = FALSE, no.margin = TRUE, label.offset = 0.001, tip.color = tipColor, edge.color = edgeVec, edge.width = edgeLwd)

dev.off()



######
## Do the same for the unconstrained tree

uncontree2 <- uncontree
uncontree2 <- drop.tip(uncontree2, 'Sphenodon_punctatus')

familyGenus <- data.frame(tip = uncontree2$tip.label, family = tax[uncontree2$tip.label, 'family'], genus = tax[uncontree2$tip.label, 'genus'], familyMono = NA, genusMono = NA)

for (i in 1:length(families)) {
	spInGroup <- tax[which(tax[, 'family'] == families[i]), 'species']
	if (length(spInGroup) > 1) {
		tmp <- extract.clade(uncontree2, getMRCA(uncontree2, spInGroup))
		if (length(spInGroup) == Ntip(tmp)) {
			familyGenus[which(familyGenus$tip %in% spInGroup), 'familyMono'] <- 'TRUE'
		} else {
			message('\t', i, ' ', families[i], ' not monophyletic.')
			familyGenus[which(familyGenus$tip %in% spInGroup), 'familyMono'] <- 'FALSE'
		}
	} else {
		familyGenus[which(familyGenus$tip %in% spInGroup), 'familyMono'] <- 'TRUE'
	}
}

# genera
for (i in 1:length(genera)) {
	spInGroup <- tax[which(tax[, 'genus'] == genera[i]), 'species']
	if (length(spInGroup) > 1) {
		tmp <- extract.clade(uncontree2, getMRCA(uncontree2, spInGroup))
		if (length(spInGroup) == Ntip(tmp)) {
			familyGenus[which(familyGenus$tip %in% spInGroup), 'genusMono'] <- 'TRUE'
		} else {
			message('\t', i, ' ', genera[i], ' not monophyletic.')
			familyGenus[which(familyGenus$tip %in% spInGroup), 'genusMono'] <- 'FALSE'
		}
	} else {
		familyGenus[which(familyGenus$tip %in% spInGroup), 'genusMono'] <- 'TRUE'
	}
}

table(familyGenus$familyMono)
table(familyGenus$genusMono)

familyGenusUncon <- familyGenus

tipColor <- rep('black', Ntip(uncontree2))

# lack of genus monophyly -> orange
# lack of family monophyly -> blue
# lack of both genus and family monophyly -> purple

tipColor[which(familyGenusUncon$genusMono == FALSE & familyGenusUncon$familyMono == TRUE)] <- 'orange'
tipColor[which(familyGenusUncon$genusMono == TRUE & familyGenusUncon$familyMono == FALSE)] <- 'blue'
tipColor[which(familyGenusUncon$genusMono == FALSE & familyGenusUncon$familyMono == FALSE)] <- 'purple'

table(tipColor)

newLabels <- uncontree2$tip.label
for (i in 1:length(newLabels)) {
	newLabels[i] <- paste0(tax[uncontree2$tip.label[i], 'family'], '_>_', newLabels[i])
}

uncontree2$tip.label <- newLabels

pdf('tier1/fulltree/bestUnconstrainedTree_monophylyColors.pdf', width = 20, height = 700)

plot.phylo(uncontree2, cex = 0.5, no.margin = TRUE, label.offset = 0.001, tip.color = tipColor)

dev.off()


# how do con and uncon compare in terms of monophyly?

table(unique(familyGenusCon[, c('family', 'familyMono')])$familyMono)
table(unique(familyGenusCon[, c('genus', 'genusMono')])$genusMono)

table(unique(familyGenusUncon[, c('family', 'familyMono')])$familyMono)
table(unique(familyGenusUncon[, c('genus', 'genusMono')])$genusMono)

# How does this compare to tier1 tree that used IQTREE-selected families?
prevIQTREE <- root(prevIQTREE, 'Sphenodon_punctatus', resolve.root = TRUE)
prevIQTREE <- ladderize(prevIQTREE)
prevIQTREE <- drop.tip(prevIQTREE, 'Sphenodon_punctatus')

familyGenus <- data.frame(tip = prevIQTREE$tip.label, family = tax[prevIQTREE$tip.label, 'family'], genus = tax[prevIQTREE$tip.label, 'genus'], familyMono = NA, genusMono = NA)

for (i in 1:length(families)) {
	spInGroup <- tax[which(tax[, 'family'] == families[i]), 'species']
	if (length(spInGroup) > 1) {
		tmp <- extract.clade(prevIQTREE, getMRCA(prevIQTREE, spInGroup))
		if (length(spInGroup) == Ntip(tmp)) {
			familyGenus[which(familyGenus$tip %in% spInGroup), 'familyMono'] <- 'TRUE'
		} else {
			message('\t', i, ' ', families[i], ' not monophyletic.')
			familyGenus[which(familyGenus$tip %in% spInGroup), 'familyMono'] <- 'FALSE'
		}
	} else {
		familyGenus[which(familyGenus$tip %in% spInGroup), 'familyMono'] <- 'TRUE'
	}
}

# genera
for (i in 1:length(genera)) {
	spInGroup <- tax[which(tax[, 'genus'] == genera[i]), 'species']
	if (length(spInGroup) > 1) {
		tmp <- extract.clade(prevIQTREE, getMRCA(prevIQTREE, spInGroup))
		if (length(spInGroup) == Ntip(tmp)) {
			familyGenus[which(familyGenus$tip %in% spInGroup), 'genusMono'] <- 'TRUE'
		} else {
			#message('\t', i, ' ', genera[i], ' not monophyletic.')
			familyGenus[which(familyGenus$tip %in% spInGroup), 'genusMono'] <- 'FALSE'
		}
	} else {
		familyGenus[which(familyGenus$tip %in% spInGroup), 'genusMono'] <- 'TRUE'
	}
}

table(unique(familyGenus[, c('family', 'familyMono')])$familyMono)
table(unique(familyGenus[, c('genus', 'genusMono')])$genusMono)





