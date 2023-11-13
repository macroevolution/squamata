# Once tier2 family-level trees have been inferred, we will have RAXML-NG recalculate all tree likelihoods. 
# This script reads in those trees, writes them as tree sets, and preps the raxml-ng HPC scripts.


if (grepl('pascal', getwd())) {
	basedir <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/'
} else {
	basedir <- '~/Dropbox/squamate_tree/phylogenetic_inference/'
}

maindir <- paste0(basedir, 'tier2')

treesDest <- paste0(basedir, 'tier2/raxmlng/trees/')

setwd(maindir)

require(ape)
require(phangorn)

# read in unique topologies from genomic bootstrapping and do some name adjustment
bs <- read.tree(paste0(basedir, 'data/all_bootstraps/IQ.ufboot.unique'))

cladeMapping <- read.csv('cladeMapping.csv', stringsAsFactors = FALSE)
fams <- colnames(cladeMapping)

# iterate through families
for (i in 1:length(fams)) {
	
	# how many unique topologies are there?
	uniqueTopo <- unique(cladeMapping[,i])
	
	for (j in 1:length(uniqueTopo)) {
		
		if (fams[i] %in% c('Colubridae', 'Gekkonidae', 'Scincidae')) {
			treeFiles <- list.files(paste0('families/', fams[i], '_uniqueTopo_', uniqueTopo[j]), pattern = '\\.treefile$', full.names = TRUE)
			
			#enforce ordering
			treeFiles <- treeFiles[order(as.numeric(gsub('_|\\.', '', stringr::str_extract(treeFiles, '_\\d\\d?\\.'))))]
			treeset <- lapply(treeFiles, read.tree)
			class(treeset) <- 'multiPhylo'
			rf <- as.numeric(RF.dist(treeset))
			
			if (all(rf == 0)) {
				treeset <- NULL	
			}

		} else {
			treeFiles <- list.files(paste0('families/', fams[i], '_uniqueTopo_', uniqueTopo[j]), pattern = '\\.runtrees$', full.names = TRUE)

			treeset <- read.tree(treeFiles)
			rf <- as.numeric(RF.dist(treeset))
			
			if (all(rf == 0)) {
				treeset <- NULL	
			}
		}

		# if treeset not NULL, write it
		if (!is.null(treeset)) {
			fn <- paste0(treesDest, fams[i], '_uniqueTopo_', j, '.trees')
			write.tree(treeset, fn)
		}	
	}
}

# write SLURM scripts

slurmHeader <- 
'#!/bin/bash
#
#SBATCH --job-name=jobName
#SBATCH --output=jobName.txt
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH -p short-28core
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=pascal.title@stonybrook.edu

cd /gpfs/scratch/ptitle/tier2/raxmlng/logs;

cmd
'


setwd(treesDest)
treesetFiles <- list.files(pattern = '\\.trees')
for (i in 1:length(treesetFiles)) {
	
	# ~/raxml-ng --evaluate --msa famConcat.aln --model ../concatenated.partitions.raxml --threads 12 --prefix evalRaxmlng_family --opt-model on --opt-branches off --tree trees.trees
	
	cmd <- paste('~/raxml-ng', '--evaluate', '--msa', paste0('../../families/', gsub('\\.trees$', '', treesetFiles[i]), '/famConcat.aln'), '--model', '../concatenated.partitions.raxml', '--threads', 16, '--opt-model', 'on', '--opt-branches', 'on', '--prefix', gsub('\\.trees$', '_raxmlEval', treesetFiles[i]), '--tree', paste0('../trees/', treesetFiles[i]))
	cmd <- paste0(cmd, ';')
		
	slurmFile <- paste0(basedir, 'tier2/raxmlng/slurm/', gsub('\\.trees$', '.slurm', treesetFiles[i]))

	tmp <- gsub('cmd', cmd, slurmHeader)
	tmp <- gsub('jobName', gsub('\\.trees$', '', treesetFiles[i]), tmp)

	write(tmp, file = slurmFile)
}
		
# execute
allSlurm <- list.files(paste0(basedir, 'tier2/raxmlng/slurm'), pattern = '\\.slurm$')
write(paste0('sbatch ', allSlurm, ';'), paste0(basedir, 'tier2/raxmlng/slurm/runTier2_eval.sh'))










