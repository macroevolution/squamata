
# Move through each gene directory, locate the alignment and build a gene tree. 

# Current method of choice: IQTREE -fast

# -----------------------------------------------------
alnDir <- '~/squam2020/fasta_March2020_alignments'
availableThreads <- parallel::detectCores()
availableThreads <- 10
iqtreePath <- '/home/ptitle@ads.iu.edu/iqtree2'

if (Sys.which(iqtreePath) == '') {
	stop("Can't find IQTREE!")
}

if (dir.exists(alnDir)) {
	setwd(alnDir)
} else {
	stop('dir not found.')
}


alnFiles <- list.files(pattern='LINS1_trimmed\\.aln$')

clusterTag <- '_onedirection.trimmedFiltered_LINS1_trimmed.aln|_onedirection.filtered_LINS1_trimmed.aln'

message(availableThreads, ' threads used.')
message(length(alnFiles), ' files found.')

numSamples <- numeric(length(alnFiles))
for (i in 1:length(alnFiles)) {
	numSamples[i] <- length(Biostrings::readDNAStringSet(alnFiles[i]))
}

# alnFiles <- alnFiles[order(numSamples, decreasing = TRUE)]
# numSamples <- sort(numSamples, decreasing = TRUE)


for (i in 1:length(alnFiles)) {
	
	message('\t', alnFiles[i])
	
	inFile <- alnFiles[i]
	nRuns <- 20

	# IQTREE2 -fast
	args <- c('-s', inFile, '-m', 'GTR+G', '-pre', gsub('\\.aln$', '', inFile), '-fast', '--runs', nRuns, '-nt', availableThreads, '--redo')

	# # IQTREE2 standard
	# nRuns <- 5
	# args <- c('-s', inFile, '-m', 'GTR+G', '-pre', gsub('\\.aln$', '', inFile), '--runs', nRuns, '-nt', availableThreads, '--redo')


	call <- system2(command = iqtreePath, args = args, stdout = TRUE, stderr = TRUE)


	# keep intermediate files?
	# keep just tree and log?
}	
		


## VERSION 2 THAT GENERATES PBS FILES FOR HPC KARST

alnDir <- '~/Dropbox/Oz_Crown_Ages/finalAlignments_June2020'
pbsDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/geneTreeInferencePBS'
qsubFile <- '~/Dropbox/squamatePhylo/2019/topologyScripts/geneTreeInferencePBS/allRuns.txt'
availableThreads <- 10

iqtreePath <- '/N/u/ptitle/Karst/iqtree2'

if (dir.exists(alnDir)) {
	setwd(alnDir)
} else {
	stop('dir not found.')
}

alnFiles <- list.files(pattern='LINS1_trimmed\\.aln$')

numSamples <- numeric(length(alnFiles))
for (i in 1:length(alnFiles)) {
	numSamples[i] <- length(Biostrings::readDNAStringSet(alnFiles[i]))
}

clusterTag <- '_onedirection.trimmedFiltered_LINS1_trimmed.aln|_onedirection.filtered_LINS1_trimmed.aln'


pbsTemplate <- 
	c(
		"#!/bin/bash",
		"#####  Constructed by HPC everywhere #####",
		"#PBS -M ptitle@iu.edu",
		"#PBS -l nodes=1:ppn=10,walltime=totalTime",
		"#PBS -l vmem=16gb",
		"#PBS -m abe",
		"#PBS -N jobName",
		"#PBS -j oe",
		"",
		"######  Module commands #####",
		"",
		"######  Job commands go below this line #####",
		"cd /N/dc2/scratch/ptitle/geneTrees;",
		"~/iqtree2 -s inFile -m GTR+G -pre prefix --runs 5 -nt 16;"
)

# totalTime <- 4:00:00:00
# jobName <- ND4_iqtree


for (i in 1:length(alnFiles)) {
	
	message('\t', alnFiles[i])
	
	inFile <- alnFiles[i]
	nRuns <- 20

	# IQTREE2 -fast
	args <- c(iqtreePath, '-s', inFile, '-m', 'GTR+G', '-pre', gsub('\\.aln$', '', inFile), '-fast', '--runs', nRuns, '-nt', availableThreads)

	# IQTREE2 standard
	nRuns <- 5
	args <- c(iqtreePath, '-s', inFile, '-m', 'GTR+G', '-pre', gsub('\\.aln$', '', inFile), '--runs', nRuns, '-nt', availableThreads)

	
	if (numSamples[i] <= 1000) {
		runTime <- "3:00:00:00"
	}
	if (numSamples[i] > 1000 & numSamples[i] <= 2000) {
		runTime <- "5:00:00:00"
	}
	if (numSamples[i] > 2000 & numSamples[i] <= 2500) {
		runTime <- "7:00:00:00"
	}
	if (numSamples[i] > 2500) {
		runTime <- "9:00:00:00"
	}
	
	gene <- gsub(clusterTag, '', inFile)
	
	
	# modify PBS template
	genePBS <- pbsTemplate
	genePBS[grep('jobName', genePBS)] <- gsub('jobName', gene, grep('jobName', genePBS, value=TRUE))
	genePBS[grep('totalTime', genePBS)] <- gsub('totalTime', runTime, grep('totalTime', genePBS, value=TRUE))
	genePBS[grep('GTR', genePBS)] <- paste0(paste(args, collapse = ' '), ';')

	write(genePBS, file = paste0(pbsDir, '/', gene, '.pbs'))

}	
		
allPBS <- list.files(pbsDir, pattern = '\\.pbs$')

write(paste0('qsub ', allPBS, ';'), file = qsubFile)





# ------------------------------------------------------------
# Version that has each tree inference run as a separate job

alnDir <- '~/Dropbox/Oz_Crown_Ages/finalAlignments_March2020'
alnDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/finalAlignments'
pbsDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/geneTreeInferencePBS'
qsubFile <- '~/Dropbox/squamatePhylo/2019/topologyScripts/geneTreeInferencePBS/allRuns.txt'
availableThreads <- 6

iqtreePath <- '/N/u/ptitle/Karst/iqtree2'

if (dir.exists(alnDir)) {
	setwd(alnDir)
} else {
	stop('dir not found.')
}

alnFiles <- list.files(pattern='LINS1_trimmed\\.aln$')

numSamples <- numeric(length(alnFiles))
for (i in 1:length(alnFiles)) {
	numSamples[i] <- length(Biostrings::readDNAStringSet(alnFiles[i]))
}

clusterTag <- '_onedirection.trimmedFiltered_LINS1_trimmed.aln|_onedirection.filtered_LINS1_trimmed.aln'


pbsTemplate <- 
	c(
		"#!/bin/bash",
		"#####  Constructed by HPC everywhere #####",
		"#PBS -M ptitle@iu.edu",
		"#PBS -l nodes=1:ppn=6,walltime=totalTime",
		"#PBS -l vmem=10gb",
		"#PBS -m abe",
		"#PBS -N jobName",
		"#PBS -j oe",
		"",
		"######  Module commands #####",
		"",
		"######  Job commands go below this line #####",
		"cd /N/dc2/scratch/ptitle/geneTrees;",
		"~/iqtree2 -s inFile -m GTR+G -pre prefix -nt 6;"
)

# totalTime <- 4:00:00:00
# jobName <- ND4_iqtree

nRuns <- 10

for (i in 1:length(alnFiles)) {
	
	message('\t', alnFiles[i])
	
	inFile <- alnFiles[i]

	# IQTREE2 -fast
	args <- c(iqtreePath, '-s', inFile, '-m', 'GTR+G', '-pre', gsub('\\.aln$', '', inFile), '-fast', '--runs', nRuns, '-nt', availableThreads)

	# IQTREE2 standard
	args <- c(iqtreePath, '-s', inFile, '-m', 'GTR+G', '-pre', gsub('\\.aln$', '_seedjj', inFile), '-seed', 'jj', '-nt', availableThreads)

	# STANDARD - for a single run
	if (numSamples[i] <= 1000) {
		runTime <- "1:00:00:00"
	}
	if (numSamples[i] > 1000 & numSamples[i] <= 2000) {
		runTime <- "1:00:00:00"
	}
	if (numSamples[i] > 2000 & numSamples[i] <= 2500) {
		runTime <- "2:00:00:00"
	}
	if (numSamples[i] > 2500) {
		runTime <- "3:00:00:00"
	}
	
	gene <- gsub(clusterTag, '', inFile)
	
	for (j in 1:nRuns) {
			
		# modify PBS template
		genePBS <- pbsTemplate
		genePBS[grep('jobName', genePBS)] <- gsub('jobName', paste0(gene, j), grep('jobName', genePBS, value=TRUE))
		genePBS[grep('totalTime', genePBS)] <- gsub('totalTime', runTime, grep('totalTime', genePBS, value=TRUE))
		genePBS[grep('GTR', genePBS)] <- gsub('jj', j, paste0(paste(args, collapse = ' '), ';'))
	
		write(genePBS, file = paste0(pbsDir, '/', gene, j, '.pbs'))
	}

}	
		
allPBS <- list.files(pbsDir, pattern = '\\.pbs$')

write(paste0('qsub ', allPBS, ';'), file = qsubFile)

# ------------------------------------------
# Examine inferred trees, identify best of 10 runs, rename as .bestTree

setwd('~/Dropbox/squamatePhylo/2019/topologyScripts/geneTreesRound1')

clusters <- list.files(pattern = 'iqtree$')
clusters <- gsub('_seed\\d\\d?', '', clusters)
clusters <- sort(unique(clusters))
clusters <- gsub('\\.iqtree$', '', clusters)

dir.create('bestDefaultRuns')
outdir <- 'bestDefaultRuns'

for (i in 1:length(clusters)) {
	
	clusterFiles <- list.files()
	clusterFiles <- clusterFiles[grep(paste0('^', clusters[i]), basename(clusterFiles))]
	logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
	treeFiles <- grep('treefile$', clusterFiles, value = TRUE) 

	defaultRuns <- c()
	for (j in 1:length(logFiles)) {
		tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
		zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
		defaultRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
	}

	sort(defaultRuns)
	bestTree <- gsub('\\.iqtree$', '.treefile', logFiles[which.max(defaultRuns)])
	
	file.copy(bestTree, paste0(outdir, '/', gsub('_seed\\d\\d?', '_bestDefault', bestTree)))
	
}


# Write PBS files for allnni round2

alnDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/finalAlignments'
pbsDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/geneTreeInferencePBS/round2'
qsubFile <- '~/Dropbox/squamatePhylo/2019/topologyScripts/geneTreeInferencePBS/round2/allRuns.txt'
availableThreads <- 6

iqtreePath <- '/N/u/ptitle/Karst/iqtree2'

if (dir.exists(alnDir)) {
	setwd(alnDir)
} else {
	stop('dir not found.')
}

alnFiles <- list.files(pattern='LINS1_trimmed\\.aln$')

numSamples <- numeric(length(alnFiles))
for (i in 1:length(alnFiles)) {
	numSamples[i] <- length(Biostrings::readDNAStringSet(alnFiles[i]))
}

clusterTag <- '_onedirection.trimmedFiltered_LINS1_trimmed.aln|_onedirection.filtered_LINS1_trimmed.aln'


pbsTemplate <- 
	c(
		"#!/bin/bash",
		"#####  Constructed by HPC everywhere #####",
		"#PBS -M ptitle@iu.edu",
		"#PBS -l nodes=1:ppn=6,walltime=totalTime",
		"#PBS -l vmem=10gb",
		"#PBS -m abe",
		"#PBS -N jobName",
		"#PBS -j oe",
		"",
		"######  Module commands #####",
		"",
		"######  Job commands go below this line #####",
		"cd /N/dc2/scratch/ptitle/geneTreesRound2;",
		"~/iqtree2 -s inFile -m GTR+G -pre prefix -nt 6 --allnni -t startingTree;"
)

# totalTime <- 4:00:00:00
# jobName <- ND4_iqtree

nRuns <- 5

for (i in 1:length(alnFiles)) {
	
	message('\t', alnFiles[i])
	
	inFile <- alnFiles[i]
	
	startTree <- gsub('\\.aln$', '_bestDefault.treefile', inFile)

	# IQTREE2 standard
	args <- c(iqtreePath, '-s', inFile, '-m', 'GTR+G', '-pre', gsub('\\.aln$', '_allnni_seedjj', inFile), '-seed', 'jj', '-nt', availableThreads, '--allnni', '-t', startTree)

	
	
	# STANDARD - for a single run
	if (numSamples[i] <= 1000) {
		runTime <- "1:00:00:00"
	}
	if (numSamples[i] > 1000 & numSamples[i] <= 2000) {
		runTime <- "1:00:00:00"
	}
	if (numSamples[i] > 2000 & numSamples[i] <= 2500) {
		runTime <- "2:00:00:00"
	}
	if (numSamples[i] > 2500) {
		runTime <- "3:00:00:00"
	}
	
	gene <- gsub(clusterTag, '', inFile)
	
	for (j in 1:nRuns) {
			
		# modify PBS template
		genePBS <- pbsTemplate
		genePBS[grep('jobName', genePBS)] <- gsub('jobName', paste0(gene, j), grep('jobName', genePBS, value=TRUE))
		genePBS[grep('totalTime', genePBS)] <- gsub('totalTime', runTime, grep('totalTime', genePBS, value=TRUE))
		genePBS[grep('GTR', genePBS)] <- gsub('jj', j, paste0(paste(args, collapse = ' '), ';'))
	
		write(genePBS, file = paste0(pbsDir, '/', gene, j, '.pbs'))
	}

}	
		
allPBS <- list.files(pbsDir, pattern = '\\.pbs$')

write(paste0('qsub ', allPBS, ';'), file = qsubFile)

# --------------------------------------------------------
# Identify the best tree out of default and allnni runs

require(stringr)
require(phangorn)

round1Dir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/geneTreesRound1'
round2Dir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/geneTreesRound2'

oldDir <- '~/tempIQTREE'

clusters <- list.files(round1Dir, pattern = 'iqtree$')
clusters <- gsub('_seed\\d\\d?', '', clusters)
clusters <- sort(unique(clusters))
clusters <- gsub('\\.iqtree$', '', clusters)

dir.create(paste0(round2Dir, '/bestTrees'))
outdir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/bestGeneTrees'

plot <- TRUE

if (plot) pdf(file = paste0(outdir, '/geneTree_IQtreeRuns.pdf'), width = 10, height = 5, onefile = TRUE)
bestLikelihoods <- numeric(length(clusters))
#prevLikelihoods <- numeric(length(clusters))
for (i in 1:length(clusters)) {
	
	locus <- gsub('_onedirection.trimmedFiltered_LINS1_trimmed|_onedirection.filtered_LINS1_trimmed', '', clusters[i])
	
	message('\t', locus)
	
	# default runs
	clusterFiles <- list.files(round1Dir, full.names = TRUE)
	clusterFiles <- clusterFiles[grep(paste0('^', clusters[i]), basename(clusterFiles))]
	logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
	treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
	labs1 <- paste0('default_', str_extract(basename(logFiles), 'seed\\d\\d?'))

	defaultRuns <- c()
	for (j in 1:length(logFiles)) {
		tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
		zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
		defaultRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
	}

	# allnni runs
	clusterFiles <- list.files(round2Dir, full.names = TRUE)
	clusterFiles <- clusterFiles[grep(paste0('^', clusters[i]), basename(clusterFiles))]	
	logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
	treeFiles2 <- grep('treefile$', clusterFiles, value = TRUE)
	labs2 <- paste0('allnni_', str_extract(basename(logFiles), 'seed\\d\\d?'))

	allnniRuns <- c()
	for (j in 1:length(logFiles)) {
		tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
		zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
		allnniRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
	}
	
	alltrees <- c(treeFiles1, treeFiles2)
	allscores <- c(defaultRuns, allnniRuns)
	names(allscores) <- c(labs1, labs2)
	sort(allscores)
		
	# which is best?
	bestFile <- alltrees[which.max(allscores)]
	
	bestLikelihoods[i] <- max(allscores)
	names(bestLikelihoods)[i] <- locus

	# # what was our best likelihood previously, for comparison
	# logFile <- list.files(oldDir, pattern = paste0('^', locus), full.name = TRUE)
	# tmp <- scan(logFile, what = 'character', sep = '\n', quiet = TRUE)
	# zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
	# prevLikelihoods[i] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)?", "\\2", zz))
	# names(prevLikelihoods)[i] <- locus
	
	# how different are the trees?
	tt <- lapply(alltrees, function(x) read.tree(x))
	rfdists <- sapply(tt, function(x) RF.dist(tt[[which.max(allscores)]], x))

	if (plot) {
		par(mfrow=c(1,2))
		plot.new()
		plot.window(xlim = c(1, length(allscores)), ylim = range(allscores))
		axis(1, at = 1:length(allscores), labels = names(allscores), las=2, cex.axis=0.7)
		axis(2, cex.axis = 0.7)
		box()
		points(1:length(allscores), allscores)
		points(which.max(allscores), allscores[which.max(allscores)], pch = 20)
		abline(v = 10.5, lty=3)
		#abline(h = prevLikelihoods[i], lwd = 0.5, lty = 4, col = gray(0.5))
		mtext('likelihood', side = 2, line = 2.2, cex = 0.75)
	
		plot.new()
		plot.window(xlim = c(1, length(allscores)), ylim = rev(range(rfdists)))
		axis(1, at = 1:length(rfdists), labels = names(allscores), las=2, cex.axis=0.7)
		axis(2, cex.axis = 0.7)
		box()
		points(1:length(rfdists), rfdists)
		points(which.max(allscores), rfdists[which.max(allscores)], pch = 20)
		abline(v = 10.5, lty=3)
		mtext('RF distance', side = 2, line = 2.2, cex = 0.75)
		mtext(locus, outer = TRUE, adj = 0.5, line = -3, font = 2)
	}
	
	file.copy(bestFile, to = paste0(outdir, '/', gsub('_allnni_seed\\d\\d?|_seed\\d\\d?', '', basename(bestFile))))
	
}

if (plot) dev.off()


# which ones were better in previous attempt?
which((bestLikelihoods < prevLikelihoods) == TRUE)








