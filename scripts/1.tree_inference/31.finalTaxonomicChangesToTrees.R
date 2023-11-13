# Final taxonomy changes to tree
# also update metadata table (add finalTipName)

# Do the same for both tier1 and tier2 trees
# Do the same for alignment

# For cases, where we are dropping one of two tips that are being merged, keep the one with more loci, regardless of name.
## If both tips are off by 1 in terms of n loci, go by total seq length.

basedir <- '~/Dropbox/Oz_Crown_Ages/'

treefileCon <- paste0(basedir, 'phylogenetic_inference/tier1/fulltree/bestTrees/fulltree_default_con_1_raxmlOpt.raxml.bestTree')
treefileUncon <- paste0(basedir, 'phylogenetic_inference/tier1/fulltree/bestTrees/fulltree_default_uncon_1_raxmlOpt.raxml.bestTree')
tier1ConstraintTreeFile <- paste0(basedir, 'phylogenetic_inference/tier1/families/fullGenomBackbone_tier1.tre')

genomdat <- read.csv(paste0(basedir, 'squamate_phylogenomics_v11_ncbi.csv'))

# read in table that allows us to interpret tip labels
metaTableFile <- paste0(basedir, 'phylogenetic_inference/fasta_March2020/metatable2.csv')
taxonTable <- read.csv(metaTableFile, stringsAsFactors=FALSE)
head(taxonTable)

# add new column: finalTipName. Taxa that are dropped will be marked here as "dropped"
taxonTable <- cbind.data.frame(taxonTable, finalTipName = taxonTable$repdbTaxon)

alnFile <- paste0(basedir, 'phylogenetic_inference/finalAlignments_June2020/renamed/concatenated.phy')


library(ape)

contree <- read.tree(treefileCon)
uncontree <- read.tree(treefileUncon)
tier1Constraint <- read.tree(tier1ConstraintTreeFile)


# read tier2 trees
tier2dir <- paste0(basedir, 'tier2/fulltrees/tier2Runs')
tier2files <- list.files(tier2dir, pattern = '\\.treefile$', full.name = TRUE)
tier2files <- tier2files[order(as.numeric(gsub('(tier2_bs)(\\d\\d?)_phase1.treefile', '\\2', basename(tier2files))))]
tier2trees <- lapply(tier2files, read.tree)

# read alignment
aln <- scan(alnFile, what = 'character', sep = '\n')

nTaxa <- as.numeric(gsub('(\\d+)(\\s)(\\d+)', '\\1', aln[1]))
nSites <- as.numeric(gsub('(\\d+)(\\s)(\\d+)', '\\3', aln[1]))

aln <- aln[-1]


# Compare data quantity to determine which tip to drop and rename.
## Drop the names that are not the name of the list element.
modlist <- list(
	'Anolis rodriguezii' = c('Anolis rodriguezii', 'Anolis rodriguezi'),
	'Aspidoscelis sexlineatus' = c('Aspidoscelis sexlineatus', 'Aspidoscelis sexlineata'),
	'Eumeces schneiderii' = c('Eumeces schneiderii', 'Eumeces schneideri'),
	'Hemidactylus triedrus' = c('Hemidactylus triedrus', 'Hemidactylus subtriedrus'),
	'Magliophis exiguus' = c('Magliophis exiguus', 'Magliophis exiguum'),
	'Mochlus sundevallii' = c('Mochlus sundevallii', 'Mochlus sundevalli'),
	'Sinonatrix percarinatus' = c('Sinonatrix percarinatus', 'Sinonatrix percarinata'),
	'Kladirostratus acutus' = c('Kladirostratus acutus', 'Psammophylax acutus'),
	'Hypsiglena ochrorhynchus' = c('Hypsiglena ochrorhynchus', 'Hypsiglena ochrorhyncha')
)

# change family classification: Psammophiidae becomes Lamprophiidae. 
nlociList <- list()
for (i in 1:length(modlist)) {
	
	nloci <- numeric(length(modlist[[i]]))
	for (j in 1:length(modlist[[i]])) {
		ind <- which(taxonTable$repdbTaxon == modlist[[i]][j])
		nloci[j] <- nrow(taxonTable[ind,])
	}
	nlociList[[i]] <- nloci	
}

for (i in 1:length(modlist)) {
	
	message('\t', i)
	
	if (abs(diff(nlociList[[i]])) < 2) {
		seqLengths <- sapply(modlist[[i]], function(x) sum(taxonTable[which(taxonTable$repdbTaxon == x), 'seqLength']))
		mostLoci <- names(which.max(seqLengths))
		leastLoci <- names(which.min(seqLengths))
	} else {
		mostLoci <- modlist[[i]][which.max(nlociList[[i]])]
		leastLoci <- modlist[[i]][which.min(nlociList[[i]])]
	}
	# taxonTable[which(taxonTable$repdbTaxon %in% gsub('_', ' ', modlist[[i]])),]
	taxonTable[which(taxonTable$repdbTaxon == mostLoci), 'finalTipName'] <- names(modlist)[i]
	taxonTable[which(taxonTable$repdbTaxon == leastLoci), 'finalTipName'] <- 'dropped'
	
	# change to tree
	contree <- drop.tip(contree, gsub('\\s+', '_', leastLoci))
	contree$tip.label[which(contree$tip.label == gsub('\\s+', '_', mostLoci))] <- gsub('\\s+', '_', names(modlist)[i])

	uncontree <- drop.tip(uncontree, gsub('\\s+', '_', leastLoci))
	uncontree$tip.label[which(uncontree$tip.label == gsub('\\s+', '_', mostLoci))] <- gsub('\\s+', '_', names(modlist)[i])
	
	# tier2
	for (j in 1:length(tier2trees)) {
		tier2trees[[j]] <- drop.tip(tier2trees[[j]], gsub('\\s+', '_', leastLoci))
		tier2trees[[j]]$tip.label[which(tier2trees[[j]]$tip.label == gsub('\\s+', '_', mostLoci))] <- gsub('\\s+', '_', names(modlist)[i])
	}
	
	# alignment
	## drop the entry for leastLoci, and rename if needed the entry for mostLoci
	ind <- grep(paste0('^', gsub('\\s+', '_', leastLoci), '\\b'), aln)
	if (length(ind) == 0) stop()
	aln <- aln[ - ind]
	ind <- grep(paste0('^', gsub('\\s+', '_', mostLoci), '\\b'), aln)
	if (length(ind) == 0) stop()
	aln[ind] <- gsub(paste0('^', gsub('\\s+', '_', mostLoci), '\\b'), gsub('\\s+', '_', names(modlist)[i]), aln[ind])
}


# single name changes
# Xenochrophis flavipunctatus renamed to Fowlea flavipunctatus
# Xenochrophis piscator renamed to Fowlea piscator
# Ptychozoon kuhli renamed to Gekko kuhli.
# Chilomeniscus stramineus renamed to Sonora straminea.
# Gehyra Cysp becomes Gehyra CYsp

singleNameChanges <- rbind(
						c('Xenochrophis flavipunctatus', 'Fowlea flavipunctatus'),
						c('Xenochrophis piscator', 'Fowlea piscator'),
						c('Ptychozoon kuhli', 'Gekko kuhli'),
						c('Chilomeniscus stramineus', 'Sonora straminea'),
						c('Gehyra Cysp', 'Gehyra CYsp'))

gsub('\\s+', '_', singleNameChanges[,1]) %in% tier1Constraint$tip.label

for (i in 1:nrow(singleNameChanges)) {
	
	message('\t', i)
	
    if (gsub('\\s+', '_', singleNameChanges[i, 1]) %in% tier1Constraint$tip.label) {
        tier1Constraint$tip.label[which(tier1Constraint$tip.label == gsub('\\s+', '_', singleNameChanges[i, 1]))] <- gsub('\\s+', '_', singleNameChanges[i, 2])
    }
    
	taxonTable[which(taxonTable$repdbTaxon == singleNameChanges[i, 1]), 'finalTipName'] <- singleNameChanges[i, 2]
	
	contree$tip.label[which(contree$tip.label == gsub('\\s+', '_', singleNameChanges[i, 1]))] <- gsub('\\s+', '_', singleNameChanges[i, 2])
	uncontree$tip.label[which(uncontree$tip.label == gsub('\\s+', '_', singleNameChanges[i, 1]))] <- gsub('\\s+', '_', singleNameChanges[i, 2])
	
	# tier2
	for (j in 1:length(tier2trees)) {
		tier2trees[[j]]$tip.label[which(tier2trees[[j]]$tip.label == gsub('\\s+', '_', singleNameChanges[i, 1]))] <- gsub('\\s+', '_', singleNameChanges[i, 2])
	}
	
	# alignment
	ind <- grep(paste0('^', gsub('\\s+', '_', singleNameChanges[i, 1]), '\\b'), aln)
	if (length(ind) == 0) stop()
	aln[ind] <- gsub(paste0('^', gsub('\\s+', '_', singleNameChanges[i, 1]), '\\b'), gsub('\\s+', '_', singleNameChanges[i, 2]), aln[ind])
}


# taxa to simply drop
taxaToDrop <- c('Correlophus belepensis', 'Eremias isfahanica', 'Asymblepharus sikimmensis', 'Lophognathus gilberti', 'Pseudoxyrhopus analabe', 'Sibon noalamina', 'Leptophis modestus')

for (i in 1:length(taxaToDrop)) {
	
	taxonTable[which(taxonTable$repdbTaxon == taxaToDrop[i]), 'finalTipName'] <- 'dropped'
	
	contree <- drop.tip(contree, gsub('\\s+', '_', taxaToDrop[i]))
	uncontree <- drop.tip(uncontree, gsub('\\s+', '_', taxaToDrop[i]))
	
	# alignment
	ind <- grep(paste0('^', gsub('\\s+', '_', taxaToDrop[i]), '\\b'), aln)
	if (length(ind) == 0) stop()
	aln <- aln[ - ind]


}
tier2trees <- lapply(tier2trees, function(x) drop.tip(x, gsub('\\s+', '_', taxaToDrop)))


# changes to higher taxonomy
taxonTable[which(taxonTable$family == 'Psammophiidae'), 'family'] <- 'Lamprophiidae'

newGenus <- apply(taxonTable, 1, function(x) {
	if (x['finalTipName'] != 'dropped') {
		strsplit(x['finalTipName'], ' ')[[1]][1]
	} else {
		x['genus']
	}
})
taxonTable$genus <- newGenus

taxonTable$finalTipName <- gsub('\\s+', '_', taxonTable$finalTipName)

# Add 'type' flag to indicate whether species have genomic data, and if so what kind. Otherwise, list as GenBank.
found <- c()
for (i in 1:nrow(genomdat)) {
	
	ind <- grep(genomdat[i, 'sample'], taxonTable$accession)
	if (length(ind) > 0) {
		taxonTable[ind, 'type'] <- genomdat[i, 'type']
		found <- c(found, i)
	}
}

taxonTable[which(is.na(taxonTable$type)), 'type'] <- 'GenBank'

# add flag for which taxa were part of genomic constraint topology
all(tier1Constraint$tip.label %in% taxonTable$finalTipName)
taxonTable$inTopoConstraint <- taxonTable$finalTipName %in% tier1Constraint$tip.label
table(unique(taxonTable[, c('finalTipName', 'inTopoConstraint')])[, 'inTopoConstraint'])

# all tips should be identical between con and uncon trees and with tier 2 tree labels.
identical(sort(contree$tip.label), sort(uncontree$tip.label))
sapply(tier2trees, function(x) identical(sort(contree$tip.label), sort(x$tip.label)))

table(grepl('\\s+', contree$tip.label))

# are tree tip names and alignment names identical?
alnNames <- sapply(strsplit(aln, ' '), function(x) x[[1]])
identical(sort(contree$tip.label), sort(alnNames))

all(contree$tip.label %in% taxonTable$finalTipName)


# write new files
write.tree(contree, paste0(basedir, 'phylogenetic_inference/tier1/fulltree/bestTrees/fulltree_default_con_1_raxmlOpt.raxml.final.tre'))
write.tree(contree, paste0(basedir, 'phylogenetic_inference/tier1/fulltree/bestTrees/fulltree_default_con_1.final.tre'))
write.tree(uncontree, paste0(basedir, 'phylogenetic_inference/tier1/fulltree/bestTrees/fulltree_default_uncon_1_raxmlOpt.raxml.final.tre'))

for (i in 1:length(tier2trees)) {
	write.tree(tier2trees[[i]], paste0(basedir, 'phylogenetic_inference/tier2/fulltrees/finaltrees/tier2_bs', i, '.final.tre'))
}

write.csv(taxonTable, paste0(basedir, 'phylogenetic_inference/fasta_March2020/metatable3.csv'), row.names = FALSE)


# write new alignment file
newAlnFile <- paste0(basedir, 'phylogenetic_inference/finalAlignments_June2020/renamed/concatenated2.phy')

write(paste0(length(aln), ' ', nSites), file = newAlnFile, append = FALSE)
write(aln, file = newAlnFile, append = TRUE)



