# Once clade-specific bootstrap trees have been inferred, the following script will move those clades into bootstrap tree specific folders, where they will then be merged and prepped for full tree final inference. Products of this script will be bootstrap-specific folders, containing a full tree constraint that contains family level constraints glued onto a genomic backbone topology. A genomic constraint is also written, which is the genomic bootstrap tree, subset to those taxa that are actually found in the full alignment. 

if (grepl('pascal', getwd())) {
	basedir <- '~/Dropbox/Oz_Crown_Ages/'
} else {
	basedir <- '~/Dropbox/squamate_tree/'
}

maindir <- paste0(basedir, 'tier2')

setwd(maindir)

require(ape)

fullAln <- paste0(basedir, 'finalAlignments_June2020/renamed/concatenated.phy')

# get all taxa in alignment
alnTaxa <- system2('pxlssq', c('-s', fullAln, '--labels'), stdout = TRUE)
head(alnTaxa)
length(alnTaxa)

# read in unique topologies from genomic bootstrapping and do some name adjustment
bs <- read.tree(paste0(basedir, 'data/all_bootstraps/IQ.ufboot.unique'))

taxonTable <- read.csv('~/Dropbox/Oz_Crown_Ages/sqCL_taxa_repDBmatching.csv', stringsAsFactors=FALSE)
# manual matching supercedes repDB auto matching, unless no manual match found
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'repDBtaxon']
head(taxonTable)
table(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch))
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'species']


for (i in 1:length(bs)) {

	bs[[i]]$node.label <- NULL

	# remove outgroups
	genomOutgroups <- c('taeGut2', 'ISIS373002', 'UMFS-10956c', 'allMis1', 'chrPic1', 'galGal5', 'hg38', 'H20145a')
	genomOutgroups %in% bs[[i]]$tip.label
	genomOutgroups <- genomOutgroups[genomOutgroups %in% bs[[i]]$tip.label]
	bs[[i]] <- root(bs[[i]], outgroup = genomOutgroups)
	bs[[i]] <- ladderize(bs[[i]])
	bs[[i]] <- drop.tip(bs[[i]], genomOutgroups)

	table(bs[[i]]$tip.label %in% taxonTable$sample) # we expect 100%
	bs[[i]]$tip.label[!bs[[i]]$tip.label %in% taxonTable$sample] # we expect zero

	newLabels <- sapply(bs[[i]]$tip.label, function(x) taxonTable[which(taxonTable$sample == x)[1], 'manualMatch'])
	table(newLabels == '' | is.na(newLabels)) # we expect none
	newLabels <- gsub('\\s+', '_', newLabels)
	names(newLabels) <- NULL
	bs[[i]]$tip.label <- newLabels

	# Drop a few taxa from genomic constraint
	bs[[i]] <- drop.tip(bs[[i]], c('Bothrocophias_hyoprora', 'Corytophanes_hernandesii', 'Lophognathus_gilberti', 'Rhineura_floridana'))
}



cladeMapping <- read.csv('cladeMapping.csv', stringsAsFactors = FALSE)

# Let's do a dry run first to search for potentially missing files
for (i in 1:length(bs)) {
	
	message('\tbs ', i)
		
	# for each clade, find the appropriate clade topology, and copy the tree over
	for (j in 1:ncol(cladeMapping)) {
		
		# identify the appropriate clade directory to pull from
		cladeTree <- list.files(paste0('families/', colnames(cladeMapping)[j], '_uniqueTopo_', cladeMapping[i, j]), pattern = '\\.treefile', full.names = TRUE)
		
		if (length(cladeTree) == 0)
		message('\t\t', colnames(cladeMapping)[j], ' -- tree not found')
	}
}


# Now run again, and move clade files into appropriate boostrap-tree-specific folders
raxmlBetter <- c()
for (i in 1:length(bs)) {
	
	dirname <- paste0('fulltrees/bs_', i)
	if (!dir.exists(dirname)) {
		dir.create(dirname)
	}
	
	# for each clade, find the appropriate clade topology, and copy the tree over
	for (j in 1:ncol(cladeMapping)) {
		
		message('\tbs ', i, ' -- ', colnames(cladeMapping)[j])
		
		# get IQTREE tree set and scores
		if (colnames(cladeMapping)[j] %in% c('Colubridae', 'Gekkonidae', 'Scincidae')) {
			
			treeFiles <- list.files(paste0('families/', colnames(cladeMapping)[j], '_uniqueTopo_', cladeMapping[i, j]), pattern = '\\.treefile', full.names = TRUE)
			logFiles <- list.files(paste0('families/', colnames(cladeMapping)[j], '_uniqueTopo_', cladeMapping[i, j]), pattern = '\\.iqtree$', full.names = TRUE)
			
			# enfore ordering
			logFiles <- logFiles[order(as.numeric(gsub('_|\\.', '', stringr::str_extract(logFiles, '_\\d\\d?\\.'))))]
			treeFiles <- treeFiles[order(as.numeric(gsub('_|\\.', '', stringr::str_extract(treeFiles, '_\\d\\d?\\.'))))]
			treeset <- lapply(treeFiles, read.tree)
			class(treeset) <- 'multiPhylo'

			# get IQTREE scores
			iqtreeScores <- c()
			for (k in 1:length(logFiles)) {
				tmp <- scan(logFiles[k], what = 'character', sep = '\n', quiet = TRUE)
				zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
				iqtreeScores[k] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
			}
			
		} else {
			
			treeFiles <- list.files(paste0('families/', colnames(cladeMapping)[j], '_uniqueTopo_', cladeMapping[i, j]), pattern = '\\.runtrees', full.names = TRUE)
			logFiles <- list.files(paste0('families/', colnames(cladeMapping)[j], '_uniqueTopo_', cladeMapping[i, j]), pattern = '\\.iqtree$', full.names = TRUE)
			treeset <- read.tree(treeFiles)
			
			# get IQTREE scores
			iqtreeScores <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
			iqtreeScores <- iqtreeScores[(grep('MULTIPLE RUNS', iqtreeScores) + 3):(grep('MAXIMUM LIKELIHOOD TREE', iqtreeScores) - 1)]
			iqtreeScores <- do.call(rbind, strsplit(iqtreeScores, '\t'))
			iqtreeScores <- as.numeric(iqtreeScores[,2])
		}
		
		# copy over top tree
		cladeDir <- paste0(dirname, '/cladeComponents')
		if (!dir.exists(cladeDir)) {
			dir.create(cladeDir)
		}

		# get raxml-ng scores
		raxmlLog <- list.files('raxmlng/logs', pattern = paste0(colnames(cladeMapping)[j], '_uniqueTopo_', cladeMapping[i, j], '_raxmlEval'), full.names = TRUE)
		
		# if none found, then that is because there was no variation in topologies, rendering this exercise unnecessary
		if (length(raxmlLog) == 1) {
			raxmlScores <- scan(raxmlLog, what = 'character', sep = '\n', quiet = TRUE)
			raxmlScores <- grep('\\#\\d\\d?, final logLikelihood', raxmlScores, value = TRUE)
			raxmlScores <- stringr::str_extract(raxmlScores, '-\\d+\\.\\d+')
			raxmlScores <- as.numeric(raxmlScores)
		} else if (length(raxmlLog) == 0) {
			raxmlScores <- NA
		} else if (length(raxmlLog) > 1) {
			stop('more than one raxml log file.')
		}
		
		if (!anyNA(raxmlScores)) {
			## how do they compare?
			# plot(iqtreeScores, raxmlScores)
			# plot(rank(iqtreeScores), rank(raxmlScores))
			# which.max(iqtreeScores)
			# which.max(raxmlScores)
			raxmlBetter <- c(raxmlBetter, which.max(iqtreeScores) != which.max(raxmlScores))
			message('\t\tsame top ranked tree: ', which.max(iqtreeScores) == which.max(raxmlScores))
		
			fn <- paste0(cladeDir, '/', colnames(cladeMapping)[j], '.tre')
			write.tree(treeset[[which.max(raxmlScores)]], fn)

		} else {
			raxmlBetter <- c(raxmlBetter, FALSE)
			message('\t\tsame top ranked tree: TRUE')
			fn <- paste0(cladeDir, '/', colnames(cladeMapping)[j], '.tre')
			write.tree(treeset[[which.max(iqtreeScores)]], fn)
		}
	}
}

table(raxmlBetter) / length(raxmlBetter)
		
		

		
		

############################
# Next step is to combine these clade trees into an overall constraint

outgroupTaxa <- c('Sphenodon_punctatus', 'Varanus_eremius', 'Strophurus_elderi', 'Kentropyx_pelviceps', 'Moloch_horridus', 'Lerista_lineopunctulata', 'Chironius_multiventris')

# read in table that allows us to interpret tip labels
## We need this to know when to properly drop outliers.
taxonTableFile <- paste0(basedir, 'fasta_March2020/metatable2.csv')
taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
taxonTable$repdbTaxon <- gsub('\\s+', '_', taxonTable$repdbTaxon)

# Corrections
taxonTable[which(taxonTable$genus %in% c('Myersophis', 'Oxyrhabdium', 'Hologerrhum', 'Cyclocorus' )), 'family'] <- 'Lamprophiidae'
taxonTable[which(taxonTable$genus == 'Xylophis'), 'family'] <- 'Pareidae'

alnTaxa <- system2('pxlssq', c('-s', fullAln, '--labels'), stdout = TRUE)
head(alnTaxa)
length(alnTaxa)

setwd(maindir)

for (i in 1:length(bs)) {
	
	message('\tbs ', i)
	
	setwd(paste0(basedir, 'tier2/fulltrees/bs_', i, '/cladeComponents'))
	cladeFiles <- list.files(pattern = '\\.tre$')
	
	# read in clade trees
	cladeList <- vector('list', length(cladeFiles))
	names(cladeList) <- gsub('\\.tre$', '', cladeFiles)

	# read in clade tree and remove outgroup taxa
	for (j in 1:length(cladeFiles)) {
				
		tmp <- read.tree(cladeFiles[j])	
	
		tax <- unique(taxonTable[which(taxonTable$family == names(cladeList)[j]), 'repdbTaxon'])
		
		toDrop <- setdiff(tmp$tip.label, tax)
		tmp <- drop.tip(tmp, toDrop)
	
		cladeList[[j]] <- tmp
	}
	
	# remove branch lengths
	for (j in 1:length(cladeList)) {
		cladeList[[j]]$edge.length <- NULL
	}
	
	# We will use this particular bootstrap tree as backbone
	contree <- bs[[i]]
	
	inGenomNotFamilies <- setdiff(contree$tip.label, unlist(lapply(cladeList, function(x) x$tip.label)))
	inGenomNotFamilies <- intersect(inGenomNotFamilies, alnTaxa)

	# subset genomic tree to one representative for each clade, in order to get the backbone
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
	
	for (j in 1:length(cladeList)) {
		fullconstraint <- bind.tree(x = fullconstraint, y = cladeList[[j]], where = which(backbone$tip.label == cladeReps[j]), position = 1)
	}
	
	# drop backbone representatives
	fullconstraint <- drop.tip(fullconstraint, setdiff(cladeReps, paste0(inGenomNotFamilies, '_backbone')))
	fullconstraint$tip.label <- gsub('_backbone', '', fullconstraint$tip.label)
	backbone$tip.label <- gsub('_backbone', '', backbone$tip.label)
	
	# plot(fullconstraint, show.tip.label=F)
	
	# did we successfully maintain the topology of the family-level backbone?
	all.equal.phylo(backbone, keep.tip(fullconstraint, backbone$tip.label))
	
	# make sure all constraint tree taxon names are found in the sequence alignment
	setdiff(fullconstraint$tip.label, alnTaxa) # we expect none
	setdiff(alnTaxa, fullconstraint$tip.label) # these are the taxa that would still be "free"	
	
	# write to disk
	write.tree(fullconstraint, paste0(basedir, 'tier2/fulltrees/bs_', i, '/bs', i, '_fullconstraint.tre'))
	

	# write genomic constraint
	contree <- keep.tip(contree, intersect(contree$tip.label, alnTaxa))
	contree <- unroot(contree)
	write.tree(contree, paste0(basedir, 'tier2/fulltrees/bs_', i, '/bs', i, '_genomConstraint.tre'))
	
}

# move constraint trees for easier moving to HPC
for (i in 1:length(bs)) {
	
	file.copy(paste0(basedir, 'tier2/fulltrees/bs_', i, '/bs', i, '_fullconstraint.tre'), paste0(basedir, 'tier2/fulltrees/allconstraints/bs', i, '_fullconstraint.tre'))
	
}


# to check, how different are these full constraint and genomic constraint trees?
## we would not expect any of them to be identical

fulltreeCon <- vector('list', length(bs))
genomicCon <- vector('list', length(bs))

for (i in 1:length(bs)) {
	setwd(paste0(basedir, 'tier2/fulltrees/bs_', i))
	fulltreeCon[[i]] <- read.tree(list.files(pattern = 'fullconstraint'))
	genomicCon[[i]] <- read.tree(list.files(pattern = 'genomConstraint'))
}

class(fulltreeCon) <- 'multiPhylo'
class(genomicCon) <- 'multiPhylo'

fulltreeRF <- phangorn::RF.dist(fulltreeCon)
genomicRF <- phangorn::RF.dist(genomicCon)

fulltreeRF <- as.matrix(fulltreeRF)
diag(fulltreeRF) <- NA
genomicRF <- as.matrix(genomicRF)
diag(genomicRF) <- NA

which(genomicRF == 0, arr.ind = T) # 4 == 16 == 17 and 27 == 29
which(fulltreeRF == 0, arr.ind = T) # 4 == 16 == 17 and 27 == 29

cladeMapping[c(4,16,17),]
cladeMapping[c(27,29),]

# SOLUTION: Simply don't run bs_16, bs_17 and bs_29 <- We won't skip them. Better reproducibility


# Generate IQTREE commands

# phase 1 will be a single run to combine all components into a fully bifurcating tree:
## iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt 12 -g bs1_fullconstraint.tre -pre tier2_bs1_phase1

# phase 2 involves providing the phase 1 tree as a starting point, but only providing the genomic tree as a constraint. This will lead to further improvements. 5 replicates would be great.
##	 If run separately, do the following:
## 		iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt 12 -g bs1_genomConstraint.tre -t tier2_bs1_phase1.treefile -pre tier2_bs1_phase2_rep1 --seed 111
##	If run sequentially, do the following:
## 		iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt 12 -g bs1_genomConstraint.tre -t tier2_bs1_phase1.treefile -pre tier2_bs1_phase2 --runs 5


# Phase 1
txtfile <- paste0(basedir, 'tier2/fulltrees/tier2_phase1.txt')


slurmHeader <- 
'#!/bin/bash
#
#SBATCH --job-name=jobName
#SBATCH --output=jobName.txt
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=167:00:00
#SBATCH -p extended-28core
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=pascal.title@stonybrook.edu

cd /gpfs/scratch/ptitle/tier2;

cmd
'

pbsHeader <- 
'#!/bin/bash
#####  Constructed by HPC everywhere #####
#PBS -M ptitle@iu.edu
#PBS -l nodes=1:ppn=12,walltime=5:00:00:00
#PBS -l vmem=150gb
#PBS -m abe
#PBS -N jobName
#PBS -j oe


cd /N/slate/ptitle/tier2;

cmd
'



# carbonate may only have 12-core nodes
# seawulf is capable of 16
nthreads <- 12

for (i in 1:length(bs)) {
			
	dir <- paste0('cd ', paste0(basedir, 'tier2/fulltrees/bs_', i), ';')
	cmd <- paste0('iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt ', nthreads, ' -g bs', i, '_fullconstraint.tre -pre tier2_bs', i, '_phase1;\n')
	
	# write(dir, file = txtfile, append = ifelse(i == 1, FALSE, TRUE))
	# write(cmd, file = txtfile, append = TRUE)
	
	pbs <- pbsHeader
	cmd <- paste0('/N/u/ptitle/Carbonate/', cmd)
	pbs <- gsub('jobName', paste0('bs', i), pbs)
	pbs <- gsub('cmd', cmd, pbs)
	
	fn <- paste0(basedir, 'tier2/fulltrees/hpcScripts/tier2_bs', i, '.pbs')
	write(pbs, fn)
}

allcmd <- list.files(paste0(basedir, 'tier2/fulltrees/hpcScripts'), pattern = '\\.pbs$')
allcmd <- allcmd[order(as.numeric(gsub('bs', '', stringr::str_extract(allcmd, 'bs\\d\\d?'))))]
allcmd <- paste0('qsub ', allcmd, ';')
write(allcmd, paste0(basedir, 'tier2/fulltrees/hpcScripts/runTier2.sh'))


# Phase 2 -- parallel
txtfile <- paste0(basedir, 'tier2/fulltrees/tier2_phase2_parallel.txt')

nthreads <- 12
nIter <- 5

for (i in 1:length(bs)) {
	
	# skip 16, 17 and 29
#	if (!i %in% c(16, 17, 29)) {
		
		dir <- paste0('cd ', paste0(basedir, 'tier2/fulltrees/bs_', i), ';')
		write(dir, file = txtfile, append = ifelse(i == 1, FALSE, TRUE))
		
		for (j in 1:nIter) {
			cmd <- paste0('iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt ', nthreads, ' -g bs', i, '_genomConstraint.tre -t tier2_bs', i, '_phase1.treefile -pre tier2_bs', i, '_phase2_rep', j,  ' --seed ', j, j, j, ';\n')
			write(cmd, file = txtfile, append = TRUE)
		}	
#	}
}

# Phase 2 -- sequential
txtfile <- paste0(basedir, 'tier2/fulltrees/tier2_phase2_sequential.txt')

nthreads <- 12
nIter <- 5

for (i in 1:length(bs)) {
	
	# skip 16, 17 and 29
#	if (!i %in% c(16, 17, 29)) {
		
		dir <- paste0('cd ', paste0(basedir, 'tier2/fulltrees/bs_', i), ';')
		
		cmd <- paste0('iqtree2 -s concatenated.phy -spp concatenated.partitions.iqtree -m GTR+G -nt ', nthreads, ' -g bs', i, '_genomConstraint.tre -t tier2_bs', i, '_phase1.treefile -pre tier2_bs', i, '_phase2 --runs ', nIter, ';\n')
		
		write(dir, file = txtfile, append = ifelse(i == 1, FALSE, TRUE))
		write(cmd, file = txtfile, append = TRUE)	
#	}
}


