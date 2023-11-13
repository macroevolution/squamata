# Tier 1: Best estimate
# Assemble separate family trees into overall squamate constraint tree
## Use genomic backbone to link up separate clades

if (grepl('pascal', getwd())) {
	basedir <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/'
} else {
	basedir <- '~/Dropbox/squamate_tree/phylogenetic_inference/'
}

bestFamilyDir <- paste0(basedir, 'tier1/fulltree/bestFamilyTrees')

fullAln <- paste0(basedir, 'finalAlignments_June2020/renamed/concatenated.phy')

# get all taxa in alignment
alnTaxa <- system2('/opt/homebrew/bin/pxlssq', c('-s', fullAln, '--labels'), stdout = TRUE)
head(alnTaxa)
length(alnTaxa)


setwd(paste0(basedir, 'tier1'))

require(ape)

outgroupTaxa <- c('Sphenodon_punctatus', 'Varanus_eremius', 'Strophurus_elderi', 'Kentropyx_pelviceps', 'Moloch_horridus', 'Lerista_lineopunctulata', 'Chironius_multiventris')

# Read in genomic constraint tree and rename tips from samples to taxon names
contree <- read.tree(paste0(basedir, 'data/all_bootstraps/all.concat_ind0.01_loci0.05_all_n5185.constraint.tre'))
contree$node.label <- NULL

# remove outgroups
genomOutgroups <- c('taeGut2', 'ISIS373002', 'UMFS-10956c', 'allMis1', 'chrPic1', 'galGal5', 'hg38', 'H20145a')
genomOutgroups %in% contree$tip.label
genomOutgroups <- genomOutgroups[genomOutgroups %in% contree$tip.label]
contree <- root(contree, outgroup = genomOutgroups)
contree <- ladderize(contree)
contree <- drop.tip(contree, genomOutgroups)


taxonTable <- read.csv('~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/sqCL_taxa_repDBmatching.csv', stringsAsFactors=FALSE)
# manual matching supercedes repDB auto matching, unless no manual match found
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'repDBtaxon']
head(taxonTable)
table(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch))
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'species']

table(contree$tip.label %in% taxonTable$sample)
contree$tip.label[!contree$tip.label %in% taxonTable$sample]

newLabels <- sapply(contree$tip.label, function(x) taxonTable[which(taxonTable$sample == x)[1], 'manualMatch'])
table(newLabels == '' | is.na(newLabels))
newLabels <- gsub('\\s+', '_', newLabels)
names(newLabels) <- NULL
contree$tip.label <- newLabels

# Drop a few taxa from genomic constraint
contree <- drop.tip(contree, c('Bothrocophias_hyoprora', 'Corytophanes_hernandesii', 'Lophognathus_gilberti', 'Rhineura_floridana'))




# read in table that allows us to interpret tip labels
## We need this to know when to properly drop outliers.
taxonTableFile <- '../fasta_March2020/metatable2.csv'
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





# CONSTRAINED

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

# remove branch lengths
for (i in 1:length(cladeList)) {
	cladeList[[i]]$edge.length <- NULL
}

inGenomNotFamilies <- setdiff(contree$tip.label, unlist(lapply(cladeList, function(x) x$tip.label)))
inGenomNotFamilies <- intersect(inGenomNotFamilies, alnTaxa)

# subset genomic tree to one representative for each clade, in order to get the backbone
## For families that don't contain any taxa that are present in the genomic constraint tree, we may be better off simply leaving those out of the constraint altogether. The alternative would be to place the family as a polytomy somewhere, but where?

cladeReps <- sapply(cladeList, function(x) sort(intersect(contree$tip.label, x$tip.label))[1])
cladesNoRep <- which(is.na(cladeReps))
cladeList <- cladeList[setdiff(1:length(cladeList), cladesNoRep)]
cladeReps <- cladeReps[setdiff(1:length(cladeReps), cladesNoRep)]
cladeReps <- c(cladeReps, inGenomNotFamilies)
backbone <- keep.tip(contree, cladeReps)
backbone$tip.label <- paste0(backbone$tip.label, '_backbone')
cladeReps <- paste0(cladeReps, '_backbone')
backbone$edge.length <- NULL
all(cladeReps %in% backbone$tip.label)

# add Sphenodon outgroup to backbone
outgroup <- list(edge=matrix(c(2, 1), 1, 2), tip.label = "Sphenodon_punctatus", Nnode=1)
class(outgroup) <- 'phylo'
backbone <- bind.tree(backbone, outgroup, where = 'root', position = 1)
backbone <- unroot(backbone)

fullconstraint <- backbone

# bind clade trees to backbone

# small test to make sure this works as expected
# testTree <- phytools::pbtree(n = 5)
# testClade <- phytools::pbtree(n = 4)
# testClade$tip.label <- c('spA','spB','spC','spD')
# testNewTree <- bind.tree(x = testTree, y = testClade, where = which(testTree$tip.label == 't3'), position = 0.2)
# par(mfrow=c(1,3))
# plot(testTree)
# plot(testClade)
# plot(testNewTree)

for (i in 1:length(cladeList)) {
	fullconstraint <- bind.tree(x = fullconstraint, y = cladeList[[i]], where = which(backbone$tip.label == cladeReps[i]), position = 1)
}

# drop backbone representatives
fullconstraint <- drop.tip(fullconstraint, setdiff(cladeReps, paste0(inGenomNotFamilies, '_backbone')))
fullconstraint$tip.label <- gsub('_backbone', '', fullconstraint$tip.label)
backbone$tip.label <- gsub('_backbone', '', backbone$tip.label)

plot(fullconstraint, show.tip.label=F)

all.equal.phylo(backbone, keep.tip(fullconstraint, backbone$tip.label))

# # make sure all constraint tree taxon names are found in the sequence alignment
setdiff(fullconstraint$tip.label, alnTaxa) # we expect none
setdiff(alnTaxa, fullconstraint$tip.label) # these are the taxa that would still be "free"


# write to disk
write.tree(fullconstraint, 'fulltree/constrained/tier1_fullconstraint.tre')


###########################################
## UNCONSTRAINED
	# # # # Since we are not using a genomic backbone, we need a backbone based purely on the supermatrix alignment. 
# # # # Here, we will identify taxa to be included in this family-level tree.
# # # # We will prioritize taxa with good locus representation.

# # # cladeReps <- character(length(clades))

# # # for (i in 1:length(clades)) {
	
	# # # tax <- unique(taxonTable[which(taxonTable$family == clades[i]), 'repdbTaxon'])
	# # # tax <- intersect(tax, cladeList[[i]]$tip.label)
	# # # zz <- taxonTable[which(taxonTable$repdbTaxon %in% tax), ]
	# # # zz <- split(zz, zz$repdbTaxon)
	# # # zz <- zz[order(sapply(zz, nrow), decreasing = TRUE)]
	# # # cladeReps[i] <- zz[[1]][1, 'repdbTaxon']
# # # }

# # # # Add Sphenodon
# # # cladeReps <- c('Sphenodon_punctatus', cladeReps)

# # # sort(sapply(cladeReps, function(x) nrow(taxonTable[which(taxonTable$repdbTaxon == x), ])))

# # # # write subset of alignment
# # # ## use phyx pxrms program to remove all but the specified taxa, then read in and write out to adjust ordering

# # # # write a file containing the taxa to keep, and then use that to subset the alignment
# # # write(cladeReps, file = paste0(basedir, 'tier1/fulltree/unconstrained/keeptaxa.txt'))
# # # args <- c('-s', fullAln, '-f', paste0(basedir, 'tier1/fulltree/unconstrained/keeptaxa.txt'), '-c', '-o', paste0(basedir, 'tier1/fulltree/unconstrained/uncon_backbone.phy'))
# # # system2('pxrms', args)
# # # qq <- readDNAStringSet(paste0(basedir, 'tier1/fulltree/unconstrained/uncon_backbone.phy'))
# # # qq <- qq[c('Sphenodon_punctatus', setdiff(names(qq), 'Sphenodon_punctatus'))]
# # # for (j in 1:length(qq)) {
	# # # write(paste0('>', names(qq)[j]), file = paste0(basedir, 'tier1/fulltree/unconstrained/uncon_backbone.aln'), append = ifelse(j == 1, FALSE, TRUE))
	# # # write(as.character(qq[[j]]), file = paste0(basedir, 'tier1/fulltree/unconstrained/uncon_backbone.aln'), append = TRUE)
# # # }
# # # file.remove(paste0(basedir, 'tier1/fulltree/unconstrained/keeptaxa.txt'))
# # # file.remove(paste0(basedir, 'tier1/fulltree/unconstrained/uncon_backbone.phy'))
# # # # system2('pxlssq', args = c('-s', paste0(basedir, 'tier1/fulltree/unconstrained/uncon_backbone.aln')))


# # # ## Family-level backbone has been inferred
# # # uncon_backbone <- read.tree(paste0(basedir, 'tier1/fulltree/unconstrained/uncon_backbone.treefile'))

# # # uncon_backbone <- drop.tip(uncon_backbone, 'Sphenodon_punctatus')
# # # cladeReps <- setdiff(uncon_backbone$tip.label, 'Sphenodon_punctatus')
# # # uncon_backbone$tip.label <- paste0(uncon_backbone$tip.label, '_backbone')
# # # cladeReps <- paste0(cladeReps, '_backbone')
# # # uncon_backbone$edge.length <- NULL
# # # all(cladeReps %in% uncon_backbone$tip.label)
# # # plot(ladderize(uncon_backbone))

# # # # add Sphenodon outgroup to backbone
# # # outgroup <- list(edge=matrix(c(2, 1), 1, 2), tip.label = "Sphenodon_punctatus", Nnode=1)
# # # class(outgroup) <- 'phylo'
# # # uncon_backbone <- bind.tree(uncon_backbone, outgroup, where = 'root', position = 1)
# # # uncon_backbone <- unroot(uncon_backbone)

# # # fullconstraint <- uncon_backbone

# # # # bind clade trees to backbone

# # # # testTree <- phytools::pbtree(n = 5)
# # # # testClade <- phytools::pbtree(n = 4)
# # # # testClade$tip.label <- c('spA','spB','spC','spD')
# # # # testNewTree <- bind.tree(x = testTree, y = testClade, where = which(testTree$tip.label == 't3'), position = 0.2)
# # # # par(mfrow=c(1,3))
# # # # plot(testTree)
# # # # plot(testClade)
# # # # plot(testNewTree)


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

# remove branch lengths
for (i in 1:length(cladeList)) {
	cladeList[[i]]$edge.length <- NULL
}

# put clade reps in same order as clade list
# cladeReps <- cladeReps[order(sapply(gsub('_backbone', '', cladeReps), function(x) which(sapply(cladeList, function(y) x %in% y$tip.label) == TRUE)))]

cladeReps <- sapply(cladeList, function(x) x$tip.label[1])
cladeReps <- paste0(cladeReps, '_backbone')

# bind together as large polytomy
fullconstraint <- phytools::starTree(cladeReps)
#fullconstraint <- multi2di(fullconstraint)
outgroup <- list(edge=matrix(c(2, 1), 1, 2), tip.label = "Sphenodon_punctatus", Nnode=1)
class(outgroup) <- 'phylo'
fullconstraint <- bind.tree(fullconstraint, outgroup, where = 'root', position = 1)
fullconstraint <- unroot(fullconstraint)

#fullconstraint$edge.length <- rep(Nedge(fullconstraint), 1e-9)

for (i in 1:length(cladeList)) {
	fullconstraint <- bind.tree(x = fullconstraint, y = cladeList[[i]], where = which(fullconstraint$tip.label == cladeReps[i]), position = 1)
}


# for (i in 1:length(cladeList)) {
	# fullconstraint <- bind.tree(x = fullconstraint, y = cladeList[[i]], where = which(uncon_backbone$tip.label == cladeReps[i]), position = 1)
# }

# drop backbone representatives
fullconstraint <- drop.tip(fullconstraint, cladeReps)
#uncon_backbone$tip.label <- gsub('_backbone', '', uncon_backbone$tip.label)

plot(ladderize(fullconstraint), show.tip.label=F)
#plot(ladderize(uncon_backbone))

for (i in 1:length(cladeList)) {
	message(all.equal.phylo(cladeList[[i]], keep.tip(fullconstraint, intersect(fullconstraint$tip.label, cladeList[[i]]$tip.label)), use.edge.length = FALSE))
}

#all.equal.phylo(uncon_backbone, keep.tip(fullconstraint, uncon_backbone$tip.label))


# # make sure all constraint tree taxon names are found in the sequence alignment
fullAln <- paste0(basedir, 'finalAlignments_June2020/renamed/concatenated.phy')

# get all taxa in alignment
alnTaxa <- system2('pxlssq', c('-s', fullAln, '--labels'), stdout = TRUE)
head(alnTaxa)
length(alnTaxa)

setdiff(fullconstraint$tip.label, alnTaxa) # we expect none
setdiff(alnTaxa, fullconstraint$tip.label) # these are the taxa that would still be "free"


# write to disk
write.tree(fullconstraint, 'fulltree/unconstrained/tier1_uncon_fullconstraint.tre')


###################################
###################################

# Compare genomic to non-constrained backbone

con_backbone <- read.tree(paste0(basedir, 'tier1/fulltree/constrained/tier1_fullconstraint.tre'))
uncon_backbone <- read.tree(paste0(basedir, 'tier1/fulltree/unconstrained/uncon_backbone.treefile'))

all(uncon_backbone$tip.label %in% con_backbone$tip.label)
setdiff(uncon_backbone$tip.label, con_backbone$tip.label)

###################################
###################################
# Write SLURM files

templateSLURM <- c(
	'#!/bin/bash',
	'#',
	'#SBATCH --job-name=jobName',
	'#SBATCH --output=outputName',
	'#SBATCH --ntasks-per-node=nCore',
	'#SBATCH --nodes=1',
	'#SBATCH --time=nHrs:00:00',
	'#SBATCH -p queue',
	'#SBATCH --mail-type=BEGIN,END',
	'#SBATCH --mail-user=pascal.title@stonybrook.edu',
	'',
	'cd workingDir;',
	'',
	'command;'
)


# For full tree inference, we will opt for:
# 	- 24 threads
# 	- no need to specify memory (although > 64 gb should be fine)
# 	- longest possible runtime: 7 days, or 168 hrs
# 	- queue: extended-28core

nIter <- 5
nCore <- 16
maxTime <- 167
queue <- 'extended-28core'

slurmDir <- paste0(basedir, 'tier1/fulltree/slurm/')

## CONSTRAINED

workingDir <- '/gpfs/scratch/ptitle/fullAln'
command <- '/gpfs/home/ptitle/iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt 16 -g tier1_fullconstraint.tre --seed iii -pre preName'

for (i in 1:nIter) {
	
	call_con <- gsub('preName', paste0('fulltree_default', '_con_', i), command)
	call_con <- gsub('iii', paste0(rep(i, 3), collapse = ''), call_con)
	
	tmpSlurm <- templateSLURM
	tmpSlurm[grep('jobName', tmpSlurm)] <- gsub('jobName', paste0('fulltree_default', '_con_', i), grep('jobName', tmpSlurm, value = TRUE))
	tmpSlurm[grep('outputName', tmpSlurm)] <- gsub('outputName', paste0('fulltree_default', '_con_', i, '.txt'), grep('outputName', tmpSlurm, value = TRUE))
	tmpSlurm[grep('nCore', tmpSlurm)] <- gsub('nCore', nCore, grep('nCore', tmpSlurm, value = TRUE))
	tmpSlurm[grep('nHrs', tmpSlurm)] <- gsub('nHrs', maxTime, grep('nHrs', tmpSlurm, value = TRUE))
	tmpSlurm[grep('queue', tmpSlurm)] <- gsub('queue', queue, grep('queue', tmpSlurm, value = TRUE))
	tmpSlurm[grep('workingDir', tmpSlurm)] <- gsub('workingDir', workingDir, grep('workingDir', tmpSlurm, value = TRUE))
	tmpSlurm[grep('command', tmpSlurm)] <- gsub('command', call_con, grep('command', tmpSlurm, value = TRUE))
	
	write(tmpSlurm, file = paste0(slurmDir, 'fulltree_con_', i, '.slurm'))
}

## UNCONSTRAINED

workingDir <- '/gpfs/scratch/ptitle/fullAln'
command <- '/gpfs/home/ptitle/iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt 16 -g tier1_uncon_fullconstraint.tre --seed iii -pre preName'

for (i in 1:nIter) {
	
	call_uncon <- gsub('preName', paste0('fulltree_default', '_uncon_', i), command)
	call_uncon <- gsub('iii', paste0(rep(i, 3), collapse = ''), call_uncon)
	
	tmpSlurm <- templateSLURM
	tmpSlurm[grep('jobName', tmpSlurm)] <- gsub('jobName', paste0('fulltree_default', '_uncon_', i), grep('jobName', tmpSlurm, value = TRUE))
	tmpSlurm[grep('outputName', tmpSlurm)] <- gsub('outputName', paste0('fulltree_default', '_uncon_', i, '.txt'), grep('outputName', tmpSlurm, value = TRUE))
	tmpSlurm[grep('nCore', tmpSlurm)] <- gsub('nCore', nCore, grep('nCore', tmpSlurm, value = TRUE))
	tmpSlurm[grep('nHrs', tmpSlurm)] <- gsub('nHrs', maxTime, grep('nHrs', tmpSlurm, value = TRUE))
	tmpSlurm[grep('queue', tmpSlurm)] <- gsub('queue', queue, grep('queue', tmpSlurm, value = TRUE))
	tmpSlurm[grep('workingDir', tmpSlurm)] <- gsub('workingDir', workingDir, grep('workingDir', tmpSlurm, value = TRUE))
	tmpSlurm[grep('command', tmpSlurm)] <- gsub('command', call_uncon, grep('command', tmpSlurm, value = TRUE))
	
	write(tmpSlurm, file = paste0(slurmDir, 'fulltree_uncon_', i, '.slurm'))
}





##### Full tree step 2: 

# once unconstrained full tree is inferred, follow up with allnni.

if (grepl('pascal', getwd())) {
	basedir <- '~/Dropbox/Oz_Crown_Ages/'
} else {
	basedir <- '~/Dropbox/squamate_tree/'
}


setwd(paste0(basedir, 'tier1'))

require(ape)

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









templateSLURM <- c(
	'#!/bin/bash',
	'#',
	'#SBATCH --job-name=jobName',
	'#SBATCH --output=outputName',
	'#SBATCH --ntasks-per-node=nCore',
	'#SBATCH --nodes=1',
	'#SBATCH --time=nHrs:00:00',
	'#SBATCH -p queue',
	'#SBATCH --mail-type=BEGIN,END',
	'#SBATCH --mail-user=pascal.title@stonybrook.edu',
	'',
	'cd workingDir;',
	'',
	'command;'
)


nIter <- 5
nCore <- 16
maxTime <- 167
queue <- 'extended-28core'

slurmDir <- paste0(basedir, 'tier1/fulltree/slurm/')


## CONSTRAINED

workingDir <- '/gpfs/scratch/ptitle/fullAln'
command <- '/gpfs/home/ptitle/iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt 16 -g genomConstraint.tre --seed iii -pre preName -t startingTree --allnni'

bestDefaultTree <- 'fulltree_default_con_2.treefile'

bestCon <- read.tree(paste0(basedir, 'tier1/fulltree/constrained/', bestDefaultTree))

all(genomTree$tip.label %in% bestCon$tip.label)
setdiff(genomTree$tip.label, bestCon$tip.label)
genomTree <- keep.tip(genomTree, intersect(genomTree$tip.label, bestCon$tip.label))
all(genomTree$tip.label %in% bestCon$tip.label)

# is topology the same?
all.equal.phylo(genomTree, keep.tip(bestCon, genomTree$tip.label), use.edge.length = FALSE)

t1 <- genomTree
t2 <- keep.tip(bestCon, genomTree$tip.label)

coTree <- phytools::cophylo(t1, t2)

pdf('cophylo.pdf', width = 20, height = 150)
plot(coTree, cex = 0.01)
dev.off()

genomConstraint <- keep.tip(bestCon, genomTree$tip.label)

genomConstraint$edge.length <- NULL
genomConstraint <- unroot(genomConstraint)

write.tree(genomConstraint, paste0(basedir, 'tier1/fulltree/constrained/genomConstraint.tre'))

for (i in 1:nIter) {
	
	call_con <- gsub('preName', paste0('fulltree_allnni', '_con_', i), command)
	call_con <- gsub('iii', paste0(rep(i, 3), collapse = ''), call_con)
	call_con <- gsub('startingTree', bestDefaultTree, call_con)
	
	tmpSlurm <- templateSLURM
	tmpSlurm[grep('jobName', tmpSlurm)] <- gsub('jobName', paste0('con_', i), grep('jobName', tmpSlurm, value = TRUE))
	tmpSlurm[grep('outputName', tmpSlurm)] <- gsub('outputName', paste0('fulltree_allnni', '_con_', i, '.txt'), grep('outputName', tmpSlurm, value = TRUE))
	tmpSlurm[grep('nCore', tmpSlurm)] <- gsub('nCore', nCore, grep('nCore', tmpSlurm, value = TRUE))
	tmpSlurm[grep('nHrs', tmpSlurm)] <- gsub('nHrs', maxTime, grep('nHrs', tmpSlurm, value = TRUE))
	tmpSlurm[grep('queue', tmpSlurm)] <- gsub('queue', queue, grep('queue', tmpSlurm, value = TRUE))
	tmpSlurm[grep('workingDir', tmpSlurm)] <- gsub('workingDir', workingDir, grep('workingDir', tmpSlurm, value = TRUE))
	tmpSlurm[grep('command', tmpSlurm)] <- gsub('command', call_con, grep('command', tmpSlurm, value = TRUE))
	
	write(tmpSlurm, file = paste0(slurmDir, 'fulltree_con_allnni_', i, '.slurm'))
}

## UNCONSTRAINED

workingDir <- '/gpfs/scratch/ptitle/fullAln'
command <- '/gpfs/home/ptitle/iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt 16 --seed iii -pre preName -t startingTree --allnni'

bestDefaultTree <- 'fulltree_default_uncon_2.treefile'

for (i in 1:nIter) {
	
	call_uncon <- gsub('preName', paste0('fulltree_allnni', '_uncon_', i), command)
	call_uncon <- gsub('iii', paste0(rep(i, 3), collapse = ''), call_uncon)
	call_uncon <- gsub('startingTree', bestDefaultTree, call_uncon)
	
	tmpSlurm <- templateSLURM
	tmpSlurm[grep('jobName', tmpSlurm)] <- gsub('jobName', paste0('uncon', i), grep('jobName', tmpSlurm, value = TRUE))
	tmpSlurm[grep('outputName', tmpSlurm)] <- gsub('outputName', paste0('fulltree_allnni', '_uncon_', i, '.txt'), grep('outputName', tmpSlurm, value = TRUE))
	tmpSlurm[grep('nCore', tmpSlurm)] <- gsub('nCore', nCore, grep('nCore', tmpSlurm, value = TRUE))
	tmpSlurm[grep('nHrs', tmpSlurm)] <- gsub('nHrs', maxTime, grep('nHrs', tmpSlurm, value = TRUE))
	tmpSlurm[grep('queue', tmpSlurm)] <- gsub('queue', queue, grep('queue', tmpSlurm, value = TRUE))
	tmpSlurm[grep('workingDir', tmpSlurm)] <- gsub('workingDir', workingDir, grep('workingDir', tmpSlurm, value = TRUE))
	tmpSlurm[grep('command', tmpSlurm)] <- gsub('command', call_uncon, grep('command', tmpSlurm, value = TRUE))
	
	write(tmpSlurm, file = paste0(slurmDir, 'fulltree_uncon_allnni_', i, '.slurm'))
}




