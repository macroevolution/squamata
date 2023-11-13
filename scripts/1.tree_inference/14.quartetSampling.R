# Run quartet sampling on each gene tree

## We will use quartet sampling (https://github.com/FePhyFoFum/quartetsampling) to identify rogue samples for each gene

setwd('~/squam2020/bestGeneTrees')

clusters <- list.files(getwd(), pattern = '\\.tre$|\\.treefile$')
clusters <- gsub('\\.tre$|\\.treefile$', '', clusters)

message('\t', length(clusters), ' trees detected.')

outdir <- '../quartetSampling'
if (!dir.exists(outdir)) {
	dir.create(outdir)
}

nCores <- parallel::detectCores()


clusterTag <- '_onedirection.trimmedFiltered_LINS1_trimmed|_onedirection.filtered_LINS1_trimmed'


# any genes to exclude? (if already run?)
# BDNF, cmos, COI, CYTB, ND2, ND4, rRNA_12S, rRNA_16S
excludeGenes <- c('^BDNF', '^cmos', '^COI', '^CYTB', '^ND2', '^ND4', '^rRNA_12S', '^rRNA_16S', '^ADNP', '^AHR', '^AKAP9', '^AMEL', '^BACH1', '^BHLHB2', '^BMP2', '^CAND1', '^CARD4', '^CILP')
excludeGenes <- c()
if (length(excludeGenes) > 0) {
	excludeGenes <- paste(excludeGenes, collapse = '|')
	clusters <- clusters[ - grep(excludeGenes, clusters)]
}

for (i in 1:length(clusters)) {

	gene <- clusters[i]
	gene <- gsub(clusterTag, '', gene)

	message('\t', gene)
		
	treeFile <- grep('\\.tre$|\\.treefile$', list.files(pattern = paste0('^', gene)), value = TRUE)
	alnFile <- grep('\\.aln$', list.files(pattern = paste0('^', gene)), value = TRUE)

	# create phylip-formatted alignment
	# ~/quartetsampling/pysrc/utils/fasta2phy.py seqs.fasta seqs.phylip
	phylipFile <- gsub('\\.aln$', '.phy', alnFile)
	args <- c('~/quartetsampling/pysrc/utils/fasta2phy.py', alnFile, phylipFile)
	system2('python3', args = args)
	
	if (length(grep('\\.phy$', list.files(pattern = gene), value = TRUE)) == 0) stop('phylip file not found.')

	# run quartetsampling
	args <- c('~/quartetsampling/pysrc/quartet_sampling.py', '--tree', treeFile, '--align', phylipFile, '--reps', '100', '--threads', nCores, '--lnlike', '2', '--result-prefix', gene, '--results-dir', outdir)
	
	system2('python3', args = args)
	
	
	# delete phylip file
	file.remove(phylipFile)
}


# VERSION FOR HPC -- run for the larger loci

# how big are the loci?
numSamples <- numeric(length(clusters))
names(numSamples) <- clusters
for (i in 1:length(clusters)) {
	gene <- clusters[i]
	gene <- gsub(clusterTag, '', gene)	
	treeFile <- grep('\\.tre$|\\.treefile$', list.files(pattern = paste0('^', gene)), value = TRUE)
	numSamples[i] <- ape::Ntip(ape::read.tree(treeFile))
}

# which are over 2000?
clusters <- clusters[which(numSamples > 2000)]
# BDNF, cmos, COI, CYTB, ND2, ND4, rRNA_12S, rRNA_16S
# BDNF already run
clusters <- clusters[ - grep('^BDNF', clusters)]


# We will set all jobs to 5 hrs

pbsTemplate <- 
	c(
		"#!/bin/bash",
		"#####  Constructed by HPC everywhere #####",
		"#PBS -M ptitle@iu.edu",
		"#PBS -l nodes=1:ppn=12,walltime=01:00:00:00",
		"#PBS -l vmem=16gb",
		"#PBS -m abe",
		"#PBS -N jobName",
		"#PBS -j oe",
		"",
		"######  Module commands #####",
		"module unload python",
		"module load python/3.6.8",
		"module unload gcc",
		"module load gcc/7.4.0",
		"module load raxmlng",
		"",
		"######  Job commands go below this line #####",
		"cd /N/dc2/scratch/ptitle/quartetSampling;",
		"program"
)

outdir <- '../qsResults'
nCores <- 12

pbsDir <- '../qsPBS'

for (i in 1:length(clusters)) {

	gene <- clusters[i]
	gene <- gsub(clusterTag, '', gene)

	message('\t', gene)
	
	treeFile <- grep('\\.tre$|\\.treefile$', list.files(pattern = gene), value = TRUE)
	alnFile <- grep('\\.aln$', list.files(pattern = gene), value = TRUE)

	# create phylip-formatted alignment
	# ~/quartetsampling/pysrc/utils/fasta2phy.py seqs.fasta seqs.phylip
	phylipFile <- gsub('\\.aln$', '.phy', alnFile)
	args <- c('~/quartetsampling/pysrc/utils/fasta2phy.py', alnFile, phylipFile)
	system2('python3', args = args)
	
	if (length(grep('\\.phy$', list.files(pattern = gene), value = TRUE)) == 0) stop('phylip file not found.')

	# run quartetsampling
	tempDir <- paste0('./QuartetSampling_', gene)
	args <- c('python3', '/N/u/ptitle/Karst/quartetsampling/pysrc/quartet_sampling.py', '--tree', treeFile, '--align', phylipFile, '--reps', '100', '--threads', nCores, '--lnlike', '2', '--result-prefix', gene, '--results-dir', outdir, '--temp-dir', tempDir)
	
	# modify PBS template
	genePBS <- pbsTemplate
	genePBS[grep('jobName', genePBS)] <- gsub('jobName', gene, grep('jobName', genePBS, value=TRUE))
	genePBS[grep('program', genePBS)] <- paste0(paste(args, collapse = ' '), ';')

	write(genePBS, file = paste0(pbsDir, '/', gene, '.pbs'))
	
	
	
}






###########################################

# Parse results

# Quartet fidelity (QF) score (for taxa):
##		proportion of total replicates across all branches that include a taxon that result in a concordant quartet topology. Ranges from 0-1. 

# Quartet concordance (QC) score (for branch):
## 		relative support among the 3 possible resolutions of 4 taxa. If most commonly sampled topology is concordant with input tree, than value is 0-1, if most commonly sampled topolgy is discordant with input tree, value is -1-0. 

# Quartet Differential (QD) score (for branch):
## 		disparity between proportions of the discordant topologies. less relevant here?

# Quartet Informativeness (QI) score (for branch):
## 		proportion of replicates where the best likelihood quartet has a likelihood that exceeds the quartet tree with second-best likelihood values. Allows one to distinguish between lack of concordance and lack of information. Values 0-1. QI of 1 = all quartets are informative, QI = 0 means there is no significant information. 


setwd('~/squam2020/fasta_March2020_alignments')
setwd('~/squam2020/bestGeneTrees')

clusters <- list.files(getwd(), pattern = '\\.tre$|\\.treefile$')
clusters <- gsub('\\.tre$|\\.treefile$', '', clusters)

clusterTag <- '_onedirection.trimmedFiltered_LINS1_trimmed|_onedirection.filtered_LINS1_trimmed'

message('\t', length(clusters), ' trees detected.')

outdir <- '../quartetSampling'

setwd(outdir)

cutoff <- 0.02

allScores <- c()

rogueList <- vector('list', length(clusters))
names(rogueList) <- gsub(clusterTag, '', clusters)
for (i in 1:length(clusters)) {

	gene <- clusters[i]
	gene <- gsub(clusterTag, '', gene)

	message('\t', gene)
	
	# get scores file
	scores <- list.files(pattern = paste0('^', gene))
	scores <- grep('node\\.scores\\.csv', scores, value = TRUE)	
	scores <- read.csv(scores, stringsAsFactors=FALSE)
	
	# how many tips?
	nTips <- length(setdiff(scores$node_label, grep('QS\\d+', scores$node_label, value = TRUE)))
	
	# here, we are primarily interested in quartet fidelity, which is equivalent to a rogue taxon test
	# We will use a QF cutoff of 0.05
	# par(mfrow=c(1,3))
	# hist(as.numeric(scores$qi), col=gray(0.95), xlab = 'quartet informativeness', main = '', breaks=50)	
	# hist(as.numeric(scores$qc), col=gray(0.95), xlab = 'quartet consistency', main = '', breaks=50)	
	# hist(as.numeric(scores$qf), col=gray(0.95), xlab = 'quartet fidelity', main = '', breaks=50)
	# abline(v= cutoff, lty=2, lwd=0.7, col='dark orange')
	
	allScores[[i]] <- scores[which(!is.na(scores$qf)), 'qf']
		
	rogueInd <- which(!is.na(scores$qf) & scores$qf < cutoff)
	if (length(rogueInd) > 0) {
		rogues <- scores[which(!is.na(scores$qf) & scores$qf < cutoff), 'node_label']
		rogueList[[i]] <- cbind.data.frame(gene = gene, nTips = nTips, qf = rogues, stringsAsFactors=FALSE)
	} else {
		message('\t\tNo rogues for ', gene, '!')
	}
}

# combine and write
rogueTaxa <- do.call(rbind, rogueList)
head(rogueTaxa)
table(rogueTaxa$gene)
sort(table(rogueTaxa$gene))
perc <- lapply(rogueList, function(x) nrow(x)/ x[1,2])
sort(round(unlist(perc), 4))
sum(table(rogueTaxa$gene))

allScores <- unlist(allScores)
# hist(allScores, breaks = 50)

write.csv(rogueTaxa, '../roguesFromGeneTrees.csv', row.names = FALSE)




















