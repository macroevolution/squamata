# family-level phylogenetic inference:
## identify taxa that will be inferred together
## generate appropriate constraint tree (and whether a constraint is necessary)
## attach appropriate taxa for rooting
## identify taxa that will not be part of family-level inference
## generate command for tree inference

## do not provide outgroup flag, but place outgroup taxa first in the alignment file

# Exclude from constraint: 
#- Bothrocophias hyoprora
#- Corytophanes hernandesii
#- Lophognathus gilberti
#- Rhineura floridiana

require(ape)
require(Biostrings)

if (grepl('pascal', getwd())) {
	basedir <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/'
} else {
	basedir <- '~/Dropbox/squamate_tree/phylogenetic_inference/'
}


familiesOutDir <- paste0(basedir, 'tier1/families/')
alnDir <- paste0(basedir, 'finalAlignments_June2020/renamed')
alnFiles <- list.files(alnDir, pattern = '\\.aln$', full.names = FALSE)

fullAln <- paste0(basedir, 'finalAlignments_June2020/renamed/concatenated.phy')

# get all taxa in alignment
alnTaxa <- system2('/opt/homebrew/bin/pxlssq', c('-s', fullAln, '--labels'), stdout = TRUE)
head(alnTaxa)
length(alnTaxa)

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


taxonTable <- read.csv('~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/sqCL_taxa_repDBmatching.csv', stringsAsFactors=FALSE)
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

# which families have 1 or 2 taxa? These will not be inferred independently
smallFamilies <- families[which(lengths(sapply(families, function(x) unique(taxonTable[taxonTable$family == x, 'repdbTaxon']))) < 3)]
excludeTaxa <- c(excludeTaxa, unique(taxonTable[which(taxonTable$family %in% smallFamilies), 'repdbTaxon']))

families <- setdiff(families, c(smallFamilies, 'Dibamidae'))

# For each family, identify the associated species, write new locus-specific alignments (in separate directories)
# Determine what the outgroup taxa should be (one from each of the 6 major clades, must have good locus representation)
# Determine what should be in the constraint tree, and whether a constraint tree is needed.
# Determine how many iterations should be run.

# identify ideal set of outgroups
# Sphenodon 
# Anguids: Varanus_eremius
# Gekko: Strophurus_elderi
# Laterata: Kentropyx_pelviceps
# Iguania: Moloch_horridus
# Skinks: Lerista_lineopunctulata
# Snakes: Chironius_multiventris
outgroupTaxa <- c('Sphenodon_punctatus', 'Varanus_eremius', 'Strophurus_elderi', 'Kentropyx_pelviceps', 'Moloch_horridus', 'Lerista_lineopunctulata', 'Chironius_multiventris')

# directory to output folders containing family-specific inference files.
outdir <- '~/Dropbox/Oz_Crown_Ages/tier1/families/'

# master text file containing IQTREE calls
allCommands <- paste0(basedir, 'tier1/family_allIQTREEcommands.txt')
allCommands_uncon <- paste0(basedir, 'tier1/family_allIQTREE_uncon_commands.txt')

pbsDir <- paste0(basedir, 'tier1/familyPBS')
pbsDirUncon <- paste0(basedir, 'tier1/familyPBS_unconstrained')
qsubFile <- paste0(basedir, 'tier1/familyPBS_submit.txt')
qsubFile_uncon <- paste0(basedir, 'tier1/familyPBS_uncon_submit.txt')

# set to only send email on aborts
pbsTemplate <-	c(
		"#!/bin/bash",
		"#####  Constructed by HPC everywhere #####",
		"#PBS -M ptitle@iu.edu",
		"#PBS -l nodes=1:ppn=12,walltime=totalTime",
		"#PBS -l vmem=30gb",
		"#PBS -m ae",
		"#PBS -N jobName",
		"#PBS -j oe",
		"",
		"######  Module commands #####",
		"",
		"######  Job commands go below this line #####",
		"cd /N/dc2/scratch/ptitle/tier1Families;",
		"cd familyPath;",
		"IQTREEcall;"
	)

# How many iterations to run?
## Smaller datasets can run faster, but larger datasets benefit from more iterations...
nIter <- 50

# We will specify 12 threads by default, but 6 for small clades (ntaxa < 50)
nCores <- 12 # is more a function of number of site patterns than number of taxa, so can remain fixed

# For smaller clades (ntaxa < 50), we can also specify 10gb memory



counts <- numeric(length(families))
names(counts) <- families
for (i in 1:length(families)) {	
	tax <- unique(taxonTable[which(taxonTable$family == families[i]), 'repdbTaxon'])
	
	# remove excluded taxa
	tax <- setdiff(tax, excludeTaxa)
	counts[i] <- length(tax)
}
sort(counts)


# toggles for output
writePBS <- TRUE
writeInputFiles <- FALSE
if (file.exists(allCommands)) file.remove(allCommands)
if (file.exists(allCommands_uncon)) file.remove(allCommands_uncon)


for (i in 1:length(families)) {
	
	tax <- unique(taxonTable[which(taxonTable$family == families[i]), 'repdbTaxon'])
	
	# are any taxa not in the concatenated alignment?
	tax <- intersect(tax, alnTaxa)
	
	# remove excluded taxa
	tax <- setdiff(tax, excludeTaxa)
	
	message('\t', families[i], ' -- ', length(tax), ' taxa')
	
	# subset constraint tree
	if (length(intersect(genomTree$tip.label, tax)) > 3) {
		con1 <- keep.tip(genomTree, intersect(genomTree$tip.label, tax))
		con1 <- unroot(con1)
		hasConstraint <- TRUE
	} else {
		con1 <- NULL
		hasConstraint <- FALSE
	}
	
	outgroups1 <- setdiff(outgroupTaxa, tax)
	
	# write files
	#- separate supermatrix for family
	#- constraint tree, if needed
	familyDir <- paste0(basedir, 'tier1/families/', families[i])

	if (writeInputFiles) {
		if (!dir.exists(familyDir)) {
			dir.create(familyDir)
		}
		
		# alignment
		## use phyx pxrms program to remove all but the specified taxa, then read in and write out to adjust ordering
		
		# write a file containing the taxa to keep, and then use that to subset the alignment
		write(c(outgroups1, tax), file = paste0(basedir, 'tier1/families/keeptaxa.txt'))
		args <- c('-s', fullAln, '-f', paste0(basedir, 'tier1/families/keeptaxa.txt'), '-c', '-o', paste0(familyDir, '/famConcat.phy'))
		system2('pxrms', args)
		qq <- readDNAStringSet(paste0(familyDir, '/famConcat.phy'))
		qq <- qq[c(outgroups1, tax)]
		for (j in 1:length(qq)) {
			write(paste0('>', names(qq)[j]), file = paste0(familyDir, '/famConcat.aln'), append = ifelse(j == 1, FALSE, TRUE))
			write(as.character(qq[[j]]), file = paste0(familyDir, '/famConcat.aln'), append = TRUE)
		}
		file.remove(paste0(basedir, 'tier1/families/keeptaxa.txt'))
		file.remove(paste0(familyDir, '/famConcat.phy'))
		# system2('pxlssq', args = c('-s', paste0(familyDir, '/famConcat.aln')))
	}

	
	# constraint tree if needed
	if (hasConstraint & writeInputFiles) {
		write.tree(con1, paste0(familyDir, '/famConstraint.tre'))
	}

	
	# Put together command
	nCores1 <- nCores
	if (length(tax) < 50) {
		nCores1 <- nCores / 2
	}
	mem1 <- 30
	if (length(tax) < 50) {
		mem1 <- 10
	}
	
	# Can set up to run families that are < 500 taxa to run sequentially (easier to manage)
	# iqtree2 -s famConcat.aln -spp concatenated.partitions.iqtree -m GTR+G -pre families[i] -g famConstraint.tre --runs niter -nt nCores
	
	args <- c('iqtree2', '-s', 'famConcat.aln', '-spp', '../concatenated.partitions.iqtree', '-m', 'GTR+G', '-nt', nCores1)
	
	if (hasConstraint) {
		args <- c(args, '-g', 'famConstraint.tre')
	}

	if (length(tax) < 500) {
		args <- c(args, '--runs', nIter, '-pre', paste0(families[i], '_default'))
	} else {
		
		args <- c(args, '--seed')
		args <- lapply(1:nIter, function(x) c(args, paste(rep(x,3), collapse = ''), '-pre', paste0(families[i], '_default_', x)))
		
	}
	
	# runtime
	# STANDARD IQTREE - for a single run
	runTime <- "2:00:00:00" # nIter sequential runs
	
	if (length(tax) > 100) { # nIter sequential runs
		runTime <- "5:00:00:00"
	}

	if (length(tax) >= 500) { # nIter independent runs
		runTime <- "4:00:00:00"
	}

	if (length(tax) >= 1000) { # nIter independent runs
		runTime <- "5:00:00:00"
	}

	# modify PBS template
	if (writePBS) {
		PBScommand <- args
		
		if (!is.list(args)) {
			PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
			tmpPBS <- pbsTemplate
			tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', families[i], grep('jobName', tmpPBS, value=TRUE))
			tmpPBS[grep('#PBS -l nodes', tmpPBS)] <- gsub('(ppn=)(\\d\\d?)', paste0('ppn=', nCores1), grep('#PBS -l nodes', tmpPBS, value = TRUE))
			tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
			tmpPBS[grep('vmem', tmpPBS)] <- gsub('(vmem=)(\\d\\d?)', paste0('vmem=', mem1), grep('vmem', tmpPBS, value=TRUE))
			tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', families[i], ';')
			tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')
	
			write(tmpPBS, file = paste0(pbsDir, '/', families[i], '.pbs'))
		} else {
			
			for (j in 1:nIter) {
				PBScommand[[j]][1] <- '/N/u/ptitle/Karst/iqtree2'
				
				tmpPBS <- pbsTemplate
				tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(families[i], j), grep('jobName', tmpPBS, value=TRUE))
				tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
				tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', families[i], ';')
				tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand[[j]], collapse = ' '), ';')
	
				write(tmpPBS, file = paste0(pbsDir, '/', families[i], '_', j, '.pbs'))
				
			}
		}

		# Write PBS file for unconstrained analysis
		## => If no constraint tree, then this is no different from constrained analysis and is not needed.
		
		if (hasConstraint) {
		
			if (!is.list(args)) {
				PBScommand <- args[!grepl('^-g$|famConstraint.tre$', args)]
				PBScommand[grep('-pre', PBScommand) + 1] <- paste0(PBScommand[grep('-pre', PBScommand) + 1], '_uncon')
	
				PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
				tmpPBS <- pbsTemplate
				tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(families[i], '_uncon'), grep('jobName', tmpPBS, value=TRUE))
				tmpPBS[grep('#PBS -l nodes', tmpPBS)] <- gsub('(ppn=)(\\d\\d?)', paste0('ppn=', nCores1), grep('#PBS -l nodes', tmpPBS, value = TRUE))
				tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
				tmpPBS[grep('vmem', tmpPBS)] <- gsub('(vmem=)(\\d\\d?)', paste0('vmem=', mem1), grep('vmem', tmpPBS, value=TRUE))
				tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', families[i], ';')
				tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')
		
				write(tmpPBS, file = paste0(pbsDirUncon, '/', families[i], '.pbs'))
			} else {
				
				for (j in 1:nIter) {
					PBScommand <- args[[j]][!grepl('^-g$|famConstraint.tre$', args[[j]])]
					PBScommand[grep('-pre', PBScommand) + 1] <- paste0(PBScommand[grep('-pre', PBScommand) + 1], '_uncon')
	
					PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
					
					tmpPBS <- pbsTemplate
					tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(families[i], '_uncon', j), grep('jobName', tmpPBS, value=TRUE))
					tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
					tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', families[i], ';')
					tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')
		
					write(tmpPBS, file = paste0(pbsDirUncon, '/', families[i], '_', j, '.pbs'))
					
				}
			}
		}
	}
	
	# add commands to master command file
	if (!is.list(args)) {
		write(paste0(paste(args, collapse = ' '), ';'), file = allCommands, append = TRUE)
	} else {
		for (j in 1:nIter) {
			write(paste0(paste(args[[j]], collapse = ' '), ';'), file = allCommands, append = TRUE)
		}
	}

	# add commands to master command file
	# unconstrained version (only if there is a constraint tree)
	if (hasConstraint) {
		if (!is.list(args)) {
			args <- args[!grepl('^-g$|famConstraint.tre$', args)]
			args[grep('-pre', args) + 1] <- paste0(args[grep('-pre', args) + 1], '_uncon')
			
			write(paste0(paste(args, collapse = ' '), ';'), file = allCommands_uncon, append = TRUE)
		} else {
			args <- lapply(args, function(x) x[!grepl('^-g$|famConstraint.tre$', x)])
			for (j in 1:length(args)) {
				args[[j]][grep('-pre', args[[j]]) + 1] <- paste0(args[[j]][grep('-pre', args[[j]]) + 1], '_uncon')
			}
			for (j in 1:nIter) {
				write(paste0(paste(args[[j]], collapse = ' '), ';'), file = allCommands_uncon, append = TRUE)
			}
		}
	}

	
	if (writeInputFiles) {
		rm(tax, con1, args, qq)
	} else {
		rm(tax, con1, args)
	}	
}

	
	
# write master submit file for PBS files
allPBS <- list.files(pbsDir, pattern = '\\.pbs$')
write(paste0('qsub ', allPBS, ';'), file = qsubFile)

allPBS <- list.files(pbsDirUncon, pattern = '\\.pbs$')
write(paste0('qsub ', allPBS, ';'), file = qsubFile_uncon)
	

##################################################################
# Confirm that all has been accounted for

includeNNI <- TRUE

for (i in 1:length(families)) {
	
	tax <- unique(taxonTable[which(taxonTable$family == families[i]), 'repdbTaxon'])
	
	# are any taxa not in the concatenated alignment?
	tax <- intersect(tax, alnTaxa)
	
	# remove excluded taxa
	tax <- setdiff(tax, excludeTaxa)
		
	# subset constraint tree
	if (length(intersect(genomTree$tip.label, tax)) > 3) {
		hasConstraint <- TRUE
	} else {
		hasConstraint <- FALSE
	}
	
	familyDir <- paste0(basedir, 'tier1/families/', families[i])

	allfiles <- list.files(familyDir)
		
	if ('famConstraint.tre' %in% allfiles & !hasConstraint) {
		message('\t', families[i], ' -- should not have constraint...')
	}
	
	if (families[i] %in% c('Gekkonidae', 'Scincidae', 'Colubridae')) {
		
		for (j in 1:50) {
			
			# constrained
			logFiles <- NA
			treeFiles <- NA	
			clusterFiles <- allfiles[grep('default', allfiles)]
			logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
			logFiles <- logFiles[!grepl('uncon', logFiles)]
			logFiles <- grep(paste0('_', j, '\\.'), logFiles, value = TRUE)
			treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
			treeFiles <- treeFiles[!grepl('uncon', treeFiles)]
			treeFiles <- grep(paste0('_', j, '\\.'), treeFiles, value = TRUE)
			
			if (length(c(logFiles, treeFiles)) != 2) {
				message('\t', families[i], ' - ', j, ' -- files missing for constrained!')
			}
			
			# unconstrained			
			logFiles <- NA
			treeFiles <- NA
			clusterFiles <- allfiles[grep('default', allfiles)]
			logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
			logFiles <- logFiles[grepl('uncon', logFiles)]
			logFiles <- grep(paste0('_', j, '_uncon\\.'), logFiles, value = TRUE)
			treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
			treeFiles <- treeFiles[grepl('uncon', treeFiles)]
			treeFiles <- grep(paste0('_', j, '_uncon\\.'), treeFiles, value = TRUE)
			
			if (length(c(logFiles, treeFiles)) != 2) {
				message('\t', families[i], ' - ', j, ' -- files missing for unconstrained!')
			}		
		}
		
		# ALLNNI
		if (includeNNI) {
			for (j in 1:10) {
				# constrained
				logFiles <- NA
				treeFiles <- NA	
				clusterFiles <- allfiles[grep('allnni', allfiles)]
				logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
				logFiles <- logFiles[!grepl('uncon', logFiles)]
				logFiles <- grep(paste0('_', j, '\\.'), logFiles, value = TRUE)
				treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
				treeFiles <- treeFiles[!grepl('uncon', treeFiles)]
				treeFiles <- grep(paste0('_', j, '\\.'), treeFiles, value = TRUE)
				
				if (length(c(logFiles, treeFiles)) != 2) {
					message('\t', families[i], ' - ', j, ' -- files missing for allnni constrained!')
				}
			
				# unconstrained			
				logFiles <- NA
				treeFiles <- NA
				clusterFiles <- allfiles[grep('allnni', allfiles)]
				logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
				logFiles <- logFiles[grepl('uncon', logFiles)]
				logFiles <- grep(paste0('_', j, '_uncon\\.'), logFiles, value = TRUE)
				treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
				treeFiles <- treeFiles[grepl('uncon', treeFiles)]
				treeFiles <- grep(paste0('_', j, '_uncon\\.'), treeFiles, value = TRUE)
				
				if (length(c(logFiles, treeFiles)) != 2) {
					message('\t', families[i], ' - ', j, ' -- files missing for allnni unconstrained!')
				}
			}
		}
		
	} else {
		
		# constrained
		logFiles <- NA
		treeFiles <- NA
		clusterFiles <- allfiles[grep('default', allfiles)]
		logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
		logFiles <- logFiles[!grepl('uncon', logFiles)]
		treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
		treeFiles <- treeFiles[!grepl('uncon', treeFiles)]
		
		if (length(c(logFiles, treeFiles)) != 2) {
			message('\t', families[i], ' -- files missing for constrained!')
		}
		
		treeTaxa <- read.tree(paste0(familyDir, '/', treeFiles))$tip.label
		if (length(intersect(treeTaxa, tax)) == 0) {
			message('\t', families[i], " -- constrained: taxa don't match.")
		}
		
		
		# unconstrained - only if there was a topological constraint
		if (hasConstraint) {
			logFiles <- NA
			treeFiles <- NA
			clusterFiles <- allfiles[grep('default', allfiles)]
			logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
			logFiles <- logFiles[grepl('uncon', logFiles)]
			treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
			treeFiles <- treeFiles[grepl('uncon', treeFiles)]	

			if (length(c(logFiles, treeFiles)) != 2) {
				message('\t', families[i], ' -- files missing for unconstrained!')
			}
			
			treeTaxa <- read.tree(paste0(familyDir, '/', treeFiles))$tip.label
			if (length(intersect(treeTaxa, tax)) == 0) {
				message('\t', families[i], " -- unconstrained: taxa don't match.")
			}

		}
		
		if (includeNNI) {
			# constrained
			logFiles <- NA
			treeFiles <- NA
			clusterFiles <- allfiles[grep('allnni', allfiles)]
			logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
			logFiles <- logFiles[!grepl('uncon', logFiles)]
			treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
			treeFiles <- treeFiles[!grepl('uncon', treeFiles)]
			
			if (length(c(logFiles, treeFiles)) != 2) {
				message('\t', families[i], ' -- files missing for allnni constrained!')
			}
			
			treeTaxa <- read.tree(paste0(familyDir, '/', treeFiles))$tip.label
			if (length(intersect(treeTaxa, tax)) == 0) {
				message('\t', families[i], " -- allnni constrained: taxa don't match.")
			}			
			
			# unconstrained - only if there was a topological constraint
			if (hasConstraint) {
				logFiles <- NA
				treeFiles <- NA
				clusterFiles <- allfiles[grep('allnni', allfiles)]
				logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
				logFiles <- logFiles[grepl('uncon', logFiles)]
				treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
				treeFiles <- treeFiles[grepl('uncon', treeFiles)]	
		
				if (length(c(logFiles, treeFiles)) != 2) {
					message('\t', families[i], ' -- files missing for allnni unconstrained!')
				}
				
				treeTaxa <- read.tree(paste0(familyDir, '/', treeFiles))$tip.label
				if (length(intersect(treeTaxa, tax)) == 0) {
					message('\t', families[i], " -- allnni unconstrained: taxa don't match.")
				}
	
			}
		}		
	}
}

##################################################################
# # Read in results and prep for allnni runs.

nIter <- 10

allCommands <- paste0(basedir, 'tier1/family_allIQTREEcommands_allnni.txt')
allCommands_uncon <- paste0(basedir, 'tier1/family_allIQTREE_uncon_commands_allnni.txt')
pbsDir <- paste0(basedir, 'tier1/familyPBS_allnni')
pbsDirUncon <- paste0(basedir, 'tier1/familyPBS_unconstrained_allnni')
qsubFile <- paste0(basedir, 'tier1/familyPBS_allnni_submit.txt')
qsubFile_uncon <- paste0(basedir, 'tier1/familyPBS_uncon_allnni_submit.txt')

if (file.exists(allCommands)) file.remove(allCommands)
if (file.exists(allCommands_uncon)) file.remove(allCommands_uncon)

writeInputFiles <- FALSE

pbsTemplate <-	c(
		"#!/bin/bash",
		"#####  Constructed by HPC everywhere #####",
		"#PBS -M ptitle@iu.edu",
		"#PBS -l nodes=1:ppn=12,walltime=totalTime",
		"#PBS -l vmem=30gb",
		"#PBS -m abe",
		"#PBS -N jobName",
		"#PBS -j oe",
		"",
		"######  Module commands #####",
		"",
		"######  Job commands go below this line #####",
		"cd /N/dc2/scratch/ptitle/tier1Families;",
		"cd familyPath;",
		"IQTREEcall;"
	)




for (i in 1:length(families)) {
	
	message('\t', families[i])
	
	familyDir <- paste0(basedir, 'tier1/families/', families[i])
	allnniFamilyDir <- paste0(basedir, 'tier1/families_allnni/', families[i])
	
	setwd(familyDir)
	
	# if not Colubridae, Gekkonidae, Scincidae, then best tree is the only .treefile
	# if one of those large families, then there is a separate treefile for each iteration and we need to identify the best one
	
	ntax <- Ntip(read.tree(list.files(pattern = '\\.treefile$')[1]))
	
	hasConstraint <- 'famConstraint.tre' %in% list.files()
	
	if (!families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
		
		bestDefaultTree <- paste0(families[i], '_default.treefile')
		
		if (hasConstraint) {
			bestDefaultTree_uncon <- paste0(families[i], '_default_uncon.treefile')
		} else {
			bestDefaultTree_uncon <- NA
		}
	} else {
		
		# CONSTRAINED
		# default runs
		clusterFiles <- list.files()
		clusterFiles <- clusterFiles[grep('default', clusterFiles)]
		logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
		logFiles <- logFiles[!grepl('uncon', logFiles)]
		treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
		treeFiles <- treeFiles[!grepl('uncon', treeFiles)]
		
		# enforce ordering
		logFiles <- logFiles[order(as.numeric(gsub('_|\\.', '', stringr::str_extract(logFiles, '_\\d\\d?\\.'))))]
		treeFiles <- treeFiles[order(as.numeric(gsub('_|\\.', '', stringr::str_extract(treeFiles, '_\\d\\d?\\.'))))]

		defaultRuns <- c()
		for (j in 1:length(logFiles)) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			defaultRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		
		bestDefaultTree	<- treeFiles[which.max(defaultRuns)]

		# UNCONSTRAINED
		# skip if there is no constraint for this family
		if (hasConstraint) {
			# default runs
			clusterFiles <- list.files()
			clusterFiles <- clusterFiles[grep('default', clusterFiles)]
			logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
			logFiles <- logFiles[grepl('uncon', logFiles)]
			treeFiles <- grep('treefile$', clusterFiles, value = TRUE)
			treeFiles <- treeFiles[grepl('uncon', treeFiles)]
			
			# enforce ordering
			logFiles <- logFiles[order(as.numeric(gsub('_|\\.', '', stringr::str_extract(logFiles, '_\\d\\d?\\.'))))]
			treeFiles <- treeFiles[order(as.numeric(gsub('_|\\.', '', stringr::str_extract(treeFiles, '_\\d\\d?\\.'))))]
	
			defaultRuns <- c()
			for (j in 1:length(logFiles)) {
				tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
				zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
				defaultRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
			}
			
			bestDefaultTree_uncon	<- treeFiles[which.max(defaultRuns)]
		} else {
			bestDefaultTree_uncon <- NA
		}		
	}	

	# Put together command
	nCores1 <- nCores
	if (ntax < 50) {
		nCores1 <- nCores / 2
	}
	mem1 <- 30
	if (ntax < 50) {
		mem1 <- 10
	}
	
	# Can set up to run families that are < 500 taxa to run sequentially (easier to manage)
	# iqtree2 -s famConcat.aln -spp concatenated.partitions.iqtree -m GTR+G -pre families[i] -g famConstraint.tre --runs niter -nt nCores --allnni -t startingTree
	
	args <- c('iqtree2', '-s', 'famConcat.aln', '-spp', '../concatenated.partitions.iqtree', '-m', 'GTR+G', '-nt', nCores1)
	
	if (hasConstraint) {
		args <- c(args, '-g', 'famConstraint.tre')
	}

	if (ntax < 500) {
		args <- c(args, '--runs', nIter, '-pre', paste0(families[i], '_allnni'), '-t', bestDefaultTree, '--allnni')
	} else {
		
		args <- c(args, '--seed')
		args <- lapply(1:nIter, function(x) c(args, paste(rep(x, 3), collapse = ''), '-pre', paste0(families[i], '_allnni_', x), '-t', bestDefaultTree, '--allnni'))
		
	}
	
	# runtime
	# STANDARD IQTREE - for a single run
	runTime <- "2:00:00:00" # nIter sequential runs
	
	if (length(tax) > 100) { # nIter sequential runs
		runTime <- "5:00:00:00"
	}

	if (length(tax) >= 500) { # nIter independent runs
		runTime <- "2:00:00:00"
	}

	if (length(tax) >= 1000) { # nIter independent runs
		runTime <- "3:00:00:00"
	}
	
	
	if (writeInputFiles) {
		# copy over input files
		if (!dir.exists(allnniFamilyDir)) {
			dir.create(allnniFamilyDir)
		}
		
		file.copy('famConcat.aln', paste0(allnniFamilyDir, '/famConcat.aln'))
		if (hasConstraint) {
			file.copy('famConstraint.tre', paste0(allnniFamilyDir, '/famConstraint.tre'))
		}
		file.copy(bestDefaultTree, paste0(allnniFamilyDir, '/', bestDefaultTree))
		if (!is.na(bestDefaultTree_uncon)) {
			file.copy(bestDefaultTree_uncon, paste0(allnniFamilyDir, '/', bestDefaultTree_uncon))
		}
	}
	

	# modify PBS template
	PBScommand <- args
	
	if (!is.list(args)) {
		PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
		tmpPBS <- pbsTemplate
		tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', families[i], grep('jobName', tmpPBS, value=TRUE))
		tmpPBS[grep('#PBS -l nodes', tmpPBS)] <- gsub('(ppn=)(\\d\\d?)', paste0('ppn=', nCores1), grep('#PBS -l nodes', tmpPBS, value = TRUE))
		tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
		tmpPBS[grep('vmem', tmpPBS)] <- gsub('(vmem=)(\\d\\d?)', paste0('vmem=', mem1), grep('vmem', tmpPBS, value=TRUE))
		tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', families[i], ';')
		tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')

		write(tmpPBS, file = paste0(pbsDir, '/', families[i], '.pbs'))
	} else {
		
		for (j in 1:nIter) {
			PBScommand[[j]][1] <- '/N/u/ptitle/Karst/iqtree2'
			
			tmpPBS <- pbsTemplate
			tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(families[i], j), grep('jobName', tmpPBS, value=TRUE))
			tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
			tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', families[i], ';')
			tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand[[j]], collapse = ' '), ';')

			write(tmpPBS, file = paste0(pbsDir, '/', families[i], '_', j, '.pbs'))
			
		}
	}
			
	# Write PBS file for unconstrained analysis
	if (hasConstraint) {
		if (!is.list(args)) {
			PBScommand <- args[!grepl('^-g$|famConstraint.tre$', args)]
			PBScommand[grep('-pre', PBScommand) + 1] <- paste0(PBScommand[grep('-pre', PBScommand) + 1], '_uncon')
			PBScommand[grep('treefile$', PBScommand)] <- bestDefaultTree_uncon
	
			PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
			tmpPBS <- pbsTemplate
			tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(families[i], '_uncon'), grep('jobName', tmpPBS, value=TRUE))
			tmpPBS[grep('#PBS -l nodes', tmpPBS)] <- gsub('(ppn=)(\\d\\d?)', paste0('ppn=', nCores1), grep('#PBS -l nodes', tmpPBS, value = TRUE))
			tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
			tmpPBS[grep('vmem', tmpPBS)] <- gsub('(vmem=)(\\d\\d?)', paste0('vmem=', mem1), grep('vmem', tmpPBS, value=TRUE))
			tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', families[i], ';')
			tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')
	
			write(tmpPBS, file = paste0(pbsDirUncon, '/', families[i], '.pbs'))
		} else {
			
			for (j in 1:nIter) {
				PBScommand[[j]] <- args[[j]][!grepl('^-g$|famConstraint.tre$', args[[j]])]
				PBScommand[[j]][grep('-pre', PBScommand[[j]]) + 1] <- paste0(PBScommand[[j]][grep('-pre', PBScommand[[j]]) + 1], '_uncon')
				PBScommand[[j]][grep('treefile$', PBScommand[[j]])] <- bestDefaultTree_uncon 
	
				PBScommand[[j]][1] <- '/N/u/ptitle/Karst/iqtree2'
				
				tmpPBS <- pbsTemplate
				tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(families[i], '_uncon', j), grep('jobName', tmpPBS, value=TRUE))
				tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
				tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', families[i], ';')
				tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand[[j]], collapse = ' '), ';')
	
				write(tmpPBS, file = paste0(pbsDirUncon, '/', families[i], '_', j, '.pbs'))
				
			}
		}
	}
	
	# add commands to master command file
	if (!is.list(args)) {
		write(paste0(paste(args, collapse = ' '), ';'), file = allCommands, append = TRUE)
	} else {
		for (j in 1:nIter) {
			write(paste0(paste(args[[j]], collapse = ' '), ';'), file = allCommands, append = TRUE)
		}
	}

	# add commands to master command file
	# unconstrained version
	
	if (hasConstraint) {
		if (!is.list(args)) {
			args <- args[!grepl('^-g$|famConstraint.tre$', args)]
			args[grep('-pre', args) + 1] <- paste0(args[grep('-pre', args) + 1], '_uncon')
			args[grep('treefile$', args)] <- bestDefaultTree_uncon
			
			write(paste0(paste(args, collapse = ' '), ';'), file = allCommands_uncon, append = TRUE)
		} else {
			args <- lapply(args, function(x) x[!grepl('^-g$|famConstraint.tre$', x)])
			for (j in 1:length(args)) {
				args[[j]][grep('-pre', args[[j]]) + 1] <- paste0(args[[j]][grep('-pre', args[[j]]) + 1], '_uncon')
				args[[j]][grep('treefile$', args[[j]])] <- bestDefaultTree_uncon 
			}
			for (j in 1:nIter) {
				write(paste0(paste(args[[j]], collapse = ' '), ';'), file = allCommands_uncon, append = TRUE)
			}
		}
	}

	rm(ntax, hasConstraint, bestDefaultTree, bestDefaultTree_uncon)
	
}


# write master submit file for PBS files
allPBS <- list.files(pbsDir, pattern = '\\.pbs$')
write(paste0('qsub ', allPBS, ';'), file = qsubFile)

allPBS <- list.files(pbsDirUncon, pattern = '\\.pbs$')
write(paste0('qsub ', allPBS, ';'), file = qsubFile_uncon)


# write separate exectute files for local runs
setwd('~/Dropbox/Oz_Crown_Ages/tier1/')

xx <- scan('~/Dropbox/Oz_Crown_Ages/tier1/family_allIQTREEcommands_allnni.txt', what='character', sep = '\n')

drops <- c('Colubridae', 'Scincidae', 'Gekkonidae', 'Agamidae', 'Dactyloidae', 'Viperidae', 'Lacertidae', 'Elapidae', 'Liolaemidae', 'Lamprophiidae', 'Chamaeleonidae', 'Gymnophthalmidae', 'Phrynosomatidae')
drops <- paste0(drops, collapse = '|')
xx <- xx[!grepl(drops, xx)]

basepath <- '~/squam2020/tier1Families/'
fams <- gsub('-pre\\s|_allnni', '', stringr::str_extract(xx, '-pre\\s.+_allnni'))

xx2 <- character(length(xx) * 2)

xx2[seq(1, length(xx2), by = 2)] <- paste0('cd ', basepath, fams, ';')
xx2[seq(2, length(xx2), by = 2)] <- xx

write(xx2, file = 'family_allIQTREEcommands_allnni_varanus.txt')



xx <- scan('~/Dropbox/Oz_Crown_Ages/tier1/family_allIQTREE_uncon_commands_allnni.txt', what='character', sep = '\n')

drops <- c('Colubridae', 'Scincidae', 'Gekkonidae', 'Agamidae', 'Dactyloidae', 'Viperidae', 'Lacertidae', 'Elapidae', 'Liolaemidae', 'Lamprophiidae', 'Chamaeleonidae', 'Gymnophthalmidae', 'Phrynosomatidae')
drops <- paste0(drops, collapse = '|')
xx <- xx[!grepl(drops, xx)]

length(xx)

basepath <- '~/squam2020/tier1Families/'
fams <- gsub('-pre\\s|_allnni', '', stringr::str_extract(xx, '-pre\\s.+_allnni'))

xx2 <- character(length(xx) * 2)

xx2[seq(1, length(xx2), by = 2)] <- paste0('cd ', basepath, fams, ';')
xx2[seq(2, length(xx2), by = 2)] <- xx

write(xx2, file = 'family_allIQTREE_uncon_commands_allnni_varanus.txt')


##################################################################
# Prep additional 5 allnni-only runs

writeFiles <- TRUE


allnniOnlyDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/allnniOnly/'
pbsDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/allnniOnly/pbs'
pbsDirUncon <- '~/Dropbox/squamatePhylo/2019/topologyScripts/allnniOnly/pbs_uncon'
qsubFile <- '~/Dropbox/squamatePhylo/2019/topologyScripts/allnniOnly/submit.sh'

nIter <- 5

pbsTemplate <-	c(
		"#!/bin/bash",
		"#####  Constructed by HPC everywhere #####",
		"#PBS -M ptitle@iu.edu",
		"#PBS -l nodes=1:ppn=12,walltime=totalTime",
		"#PBS -l vmem=30gb",
		"#PBS -m abe",
		"#PBS -N jobName",
		"#PBS -j oe",
		"",
		"######  Module commands #####",
		"",
		"######  Job commands go below this line #####",
		"cd /N/slate/ptitle/allnniOnly;",
		"cd familyPath;",
		"IQTREEcall;"
	)


for (i in 1:length(families)) {
	
	tax <- unique(taxonTable[which(taxonTable$family == families[i]), 'repdbTaxon'])
	
	# are any taxa not in the concatenated alignment?
	tax <- intersect(tax, alnTaxa)
	
	# remove excluded taxa
	tax <- setdiff(tax, excludeTaxa)
		
	# subset constraint tree
	if (length(intersect(genomTree$tip.label, tax)) > 3) {
		hasConstraint <- TRUE
	} else {
		hasConstraint <- FALSE
	}
		
	message('\t', families[i])
	
	familyDir <- paste0(basedir, 'tier1/families/', families[i])
	
	setwd(familyDir)
		
	destDir <- paste0(allnniOnlyDir, families[i], '/')
	if (!dir.exists(destDir)) {
		dir.create(destDir)
	}
	
	if (writeFiles) {
		file.copy('famConcat.aln', to = destDir)
		if (hasConstraint) {
			file.copy('famConstraint.tre', to = destDir)
		}
	}	
	
	# Put together command
	nCores1 <- nCores
	if (length(tax) < 50) {
		nCores1 <- nCores / 2
	}
	mem1 <- 30
	if (length(tax) < 50) {
		mem1 <- 10
	}
	
	# Can set up to run families that are < 500 taxa to run sequentially (easier to manage)
	# iqtree2 -s famConcat.aln -spp concatenated.partitions.iqtree -m GTR+G -pre families[i] -g famConstraint.tre --runs niter -nt nCores --allnni -t startingTree
	
	args <- c('iqtree2', '-s', 'famConcat.aln', '-spp', '../concatenated.partitions.iqtree', '-m', 'GTR+G', '-nt', nCores1)
	
	if (hasConstraint) {
		args <- c(args, '-g', 'famConstraint.tre')
	}
	
	# add allnni arguments
	args <- c(args, '--allnni')

	if (length(tax) < 300) {
		args <- c(args, '--runs', nIter, '-pre', paste0(families[i], '_allniOnly'))
	} else {
		
		args <- c(args, '--seed')
		args <- lapply(1:nIter, function(x) c(args, paste(rep(x, 3), collapse = ''), '-pre', paste0(families[i], '_allnniOnly_', x)))
		
	}
	
	# runtime
	# STANDARD IQTREE - for a single run
	runTime <- "2:00:00:00" # nIter sequential runs
	
	if (length(tax) > 100) { # nIter sequential runs
		runTime <- "5:00:00:00"
	}

	if (length(tax) >= 300) { # nIter independent runs
		runTime <- "3:00:00:00"
	}

	if (length(tax) >= 1000) { # nIter independent runs
		runTime <- "4:00:00:00"
	}

	# modify PBS template
	PBScommand <- args
	
	if (!is.list(args)) {
		PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
		tmpPBS <- pbsTemplate
		tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', families[i], grep('jobName', tmpPBS, value=TRUE))
		tmpPBS[grep('#PBS -l nodes', tmpPBS)] <- gsub('(ppn=)(\\d\\d?)', paste0('ppn=', nCores1), grep('#PBS -l nodes', tmpPBS, value = TRUE))
		tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
		tmpPBS[grep('vmem', tmpPBS)] <- gsub('(vmem=)(\\d\\d?)', paste0('vmem=', mem1), grep('vmem', tmpPBS, value=TRUE))
		tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./', families[i], ';')
		tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')

		write(tmpPBS, file = paste0(pbsDir, '/', families[i], '_con.pbs'))
	} else {
		
		for (j in 1:nIter) {
			PBScommand[[j]][1] <- '/N/u/ptitle/Karst/iqtree2'
			
			tmpPBS <- pbsTemplate
			tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(families[i], j), grep('jobName', tmpPBS, value=TRUE))
			tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
			tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./', families[i], ';')
			tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand[[j]], collapse = ' '), ';')

			write(tmpPBS, file = paste0(pbsDir, '/', families[i], '_', j, '_con.pbs'))
			
		}
	}
			
	# Write PBS file for unconstrained analysis
	if (hasConstraint) {
		if (!is.list(args)) {
			PBScommand <- args[!grepl('^-g$|famConstraint.tre$', args)]
			PBScommand[grep('-pre', PBScommand) + 1] <- paste0(PBScommand[grep('-pre', PBScommand) + 1], '_uncon')
	
			PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
			tmpPBS <- pbsTemplate
			tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(families[i], '_uncon'), grep('jobName', tmpPBS, value=TRUE))
			tmpPBS[grep('#PBS -l nodes', tmpPBS)] <- gsub('(ppn=)(\\d\\d?)', paste0('ppn=', nCores1), grep('#PBS -l nodes', tmpPBS, value = TRUE))
			tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
			tmpPBS[grep('vmem', tmpPBS)] <- gsub('(vmem=)(\\d\\d?)', paste0('vmem=', mem1), grep('vmem', tmpPBS, value=TRUE))
			tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./', families[i], ';')
			tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')
	
			write(tmpPBS, file = paste0(pbsDir, '/', families[i], '_uncon.pbs'))
		} else {
			
			for (j in 1:nIter) {
				PBScommand[[j]] <- args[[j]][!grepl('^-g$|famConstraint.tre$', args[[j]])]
				PBScommand[[j]][grep('-pre', PBScommand[[j]]) + 1] <- paste0(PBScommand[[j]][grep('-pre', PBScommand[[j]]) + 1], '_uncon')
	
				PBScommand[[j]][1] <- '/N/u/ptitle/Karst/iqtree2'
				
				tmpPBS <- pbsTemplate
				tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(families[i], '_uncon', j), grep('jobName', tmpPBS, value=TRUE))
				tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
				tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./', families[i], ';')
				tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand[[j]], collapse = ' '), ';')
	
				write(tmpPBS, file = paste0(pbsDir, '/', families[i], '_', j, '_uncon.pbs'))
				
			}
		}
	}	
}

allPBS <- list.files(pbsDir, pattern = '\\.pbs$')
write(paste0('qsub ', allPBS, ';'), file = qsubFile)



## Move them into the main directory
loc1 <- allnniOnlyDir
loc2 <- paste0(basedir, 'tier1/families/')

for (i in 1:length(families)) {
	
	loc1a <- paste0(loc1, families[i])
	loc2a <- paste0(loc2, families[i])
	
	f <- list.files(loc1a, full.names = TRUE)
	f <- f[grep('allnn?iOnly', basename(f))]
	
	if (!any(grepl('allnn?iOnly', list.files(loc2a)))) {
		file.copy(f, loc2a)
	}
}






##################################################################
# Prep additional thorough runs

writeFiles <- TRUE

thoroughClades <- c('Colubridae', 'Gekkonidae', 'Scincidae')
thoroughClades <- families

thoroughDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/thoroughRuns/'
pbsDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/thoroughRuns/pbs'
pbsDirUncon <- '~/Dropbox/squamatePhylo/2019/topologyScripts/thoroughRuns/pbs_uncon'
qsubFile <- '~/Dropbox/squamatePhylo/2019/topologyScripts/thoroughRuns/submit.sh'

nIter <- 5

pbsTemplate <-	c(
		"#!/bin/bash",
		"#####  Constructed by HPC everywhere #####",
		"#PBS -M ptitle@iu.edu",
		"#PBS -l nodes=1:ppn=12,walltime=totalTime",
		"#PBS -l vmem=30gb",
		"#PBS -m abe",
		"#PBS -N jobName",
		"#PBS -j oe",
		"",
		"######  Module commands #####",
		"",
		"######  Job commands go below this line #####",
		"cd /N/dc2/scratch/ptitle/tier1Families;",
		"cd familyPath;",
		"IQTREEcall;"
	)


for (i in 1:length(thoroughClades)) {
	
	tax <- unique(taxonTable[which(taxonTable$family == thoroughClades[i]), 'repdbTaxon'])
	
	# are any taxa not in the concatenated alignment?
	tax <- intersect(tax, alnTaxa)
	
	# remove excluded taxa
	tax <- setdiff(tax, excludeTaxa)
		
	# subset constraint tree
	if (length(intersect(genomTree$tip.label, tax)) > 3) {
		hasConstraint <- TRUE
	} else {
		hasConstraint <- FALSE
	}
		
	message('\t', thoroughClades[i])
	
	familyDir <- paste0(basedir, 'tier1/families/', thoroughClades[i])
	
	setwd(familyDir)
	
	# CONSTRAINED
	# default runs
	clusterFiles <- list.files()
	clusterFiles <- clusterFiles[grep('default', clusterFiles)]
	logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
	logFiles <- logFiles[!grepl('uncon', logFiles)]
	treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
	treeFiles1 <- treeFiles1[!grepl('uncon', treeFiles1)]
	iterTreeFiles <- grep('runtrees$', clusterFiles, value = TRUE)
	iterTreeFiles <- iterTreeFiles[!grepl('uncon', iterTreeFiles)]
	
	# parallel runs
	if (thoroughClades[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
		
		# enforce ordering
		logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
		treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]
		
		defaultRuns <- numeric(50)
		for (j in 1:50) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			defaultRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		iterTreeFiles <- treeFiles1
		
	# sequential runs	
	} else {
		if (length(logFiles) > 1) stop()
		labs1 <- paste0('default_constrained_', thoroughClades[i])
		tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
		zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
		zz <- do.call(rbind, strsplit(zz, '\t'))
		defaultRuns <- as.numeric(zz[,2])
	}	
	
	if (hasConstraint) {
		# UNCONSTRAINED
		# default runs
		clusterFiles <- list.files()
		clusterFiles <- clusterFiles[grep('default', clusterFiles)]
		logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
		logFiles <- logFiles[grepl('uncon', logFiles)]
		treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
		treeFiles1 <- treeFiles1[grepl('uncon', treeFiles1)]
		iterTreeFiles_uncon <- grep('runtrees$', clusterFiles, value = TRUE)
		iterTreeFiles_uncon <- iterTreeFiles_uncon[grepl('uncon', iterTreeFiles_uncon)]
		
		# parallel runs
		if (thoroughClades[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {

		# enforce ordering
		logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
		treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]

		defaultRuns_uncon <- numeric(50)
		for (j in 1:50) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			defaultRuns_uncon[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		iterTreeFiles_uncon <- treeFiles1			
			
		# sequential runs	
		} else {
			if (length(logFiles) > 1) stop()
			labs1 <- paste0('default_constrained_', thoroughClades[i])
			tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
			zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
			zz <- do.call(rbind, strsplit(zz, '\t'))
			defaultRuns_uncon <- as.numeric(zz[,2])	
		}
	} else {
		defaultRuns_uncon <- NA
	}
		
	# ALLNNI runs
	
	# CONSTRAINED
	clusterFiles <- list.files()
	clusterFiles <- clusterFiles[grep('allnni', clusterFiles)]
	logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
	logFiles <- logFiles[!grepl('uncon', logFiles)]
	treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
	treeFiles1 <- treeFiles1[!grepl('uncon', treeFiles1)]
	iterTreeFiles_allnni <- grep('runtrees$', clusterFiles, value = TRUE)
	iterTreeFiles_allnni <- iterTreeFiles_allnni[!grepl('uncon', iterTreeFiles_allnni)]

	# parallel runs
	if (thoroughClades[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {

		# enforce ordering
		logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
		treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]
		
		allnniRuns <- numeric(10)
		for (j in 1:10) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			allnniRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		iterTreeFiles_allnni <- treeFiles1
		
	# sequential runs	
	} else {
		if (length(logFiles) > 1) stop()
		labs1 <- paste0('default_constrained_', thoroughClades[i])
		tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
		zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
		zz <- do.call(rbind, strsplit(zz, '\t'))
		allnniRuns <- as.numeric(zz[,2])
	}	
	
	if (hasConstraint) {
		# UNCONSTRAINED
		# allnni runs
		clusterFiles <- list.files()
		clusterFiles <- clusterFiles[grep('allnni', clusterFiles)]
		logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
		logFiles <- logFiles[grepl('uncon', logFiles)]
		treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
		treeFiles1 <- treeFiles1[grepl('uncon', treeFiles1)]
		iterTreeFiles_allnni_uncon <- grep('runtrees$', clusterFiles, value = TRUE)
		iterTreeFiles_allnni_uncon <- iterTreeFiles_allnni_uncon[grepl('uncon', iterTreeFiles_allnni_uncon)]
		
		# parallel runs
		if (thoroughClades[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {

		# enforce ordering
		logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
		treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]

		allnniRuns_uncon <- numeric(10)
		for (j in 1:10) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			allnniRuns_uncon[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		iterTreeFiles_allnni_uncon <- treeFiles1
			
		# sequential runs	
		} else {
			if (length(logFiles) > 1) stop()
			labs1 <- paste0('default_constrained_', thoroughClades[i])
			tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
			zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
			zz <- do.call(rbind, strsplit(zz, '\t'))
			allnniRuns_uncon <- as.numeric(zz[,2])	
		}
	} else {
		allnniRuns_uncon <- NA
	}
	
	destDir <- paste0(thoroughDir, thoroughClades[i], '/')
	if (!dir.exists(destDir)) {
		dir.create(destDir)
	}
	
	if (writeFiles) {
		file.copy('famConcat.aln', to = destDir)
		if (hasConstraint) {
			file.copy('famConstraint.tre', to = destDir)
		}
	}
	
	# how different are the trees?
	# constrained
	loglik_con <- c(defaultRuns, allnniRuns)
	if (thoroughClades[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
		iterTrees <- c(as.list(iterTreeFiles), as.list(iterTreeFiles_allnni))
		if (writeFiles) {
			bestInd <- which(loglik_con == max(loglik_con))
			file.copy(iterTrees[[bestInd]], to = paste0(destDir, thoroughClades[i], '_con.treefile'))
		}

	} else {
		if (writeFiles) {
			bestInd <- which(loglik_con == max(loglik_con))
			iterTrees <- c(read.tree(iterTreeFiles), read.tree(iterTreeFiles_allnni))
			write.tree(iterTrees[[bestInd]], file = paste0(destDir, thoroughClades[i], '_con.treefile'))
		}
	}
	

	bestDefaultTree <- paste0(thoroughClades[i], '_con.treefile')
	
	# unconstrained
	if (hasConstraint) {
		loglik_uncon <- c(defaultRuns_uncon, allnniRuns_uncon)
		if (thoroughClades[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
			iterTrees_uncon <- c(as.list(iterTreeFiles_uncon), as.list(iterTreeFiles_allnni_uncon))
			if (writeFiles) {
				bestInd <- which(loglik_uncon == max(loglik_uncon))
				file.copy(iterTrees_uncon[[bestInd]], to = paste0(destDir, thoroughClades[i], '_uncon.treefile'))
			}
	
		} else {
			if (writeFiles) {
				bestInd <- which(loglik_uncon == max(loglik_uncon))
				iterTrees_uncon <- c(read.tree(iterTreeFiles_uncon), read.tree(iterTreeFiles_allnni_uncon))
				write.tree(iterTrees_uncon[[bestInd]], file = paste0(destDir, thoroughClades[i], '_uncon.treefile'))
			}
		}
				
		bestDefaultTree_uncon <- paste0(thoroughClades[i], '_uncon.treefile')
	} 
	
	# Put together command
	nCores1 <- nCores
	if (length(tax) < 50) {
		nCores1 <- nCores / 2
	}
	mem1 <- 30
	if (length(tax) < 50) {
		mem1 <- 10
	}
	
	# Can set up to run families that are < 500 taxa to run sequentially (easier to manage)
	# iqtree2 -s famConcat.aln -spp concatenated.partitions.iqtree -m GTR+G -pre families[i] -g famConstraint.tre --runs niter -nt nCores --allnni -t startingTree
	
	args <- c('iqtree2', '-s', 'famConcat.aln', '-spp', '../concatenated.partitions.iqtree', '-m', 'GTR+G', '-nt', nCores1)
	
	if (hasConstraint) {
		args <- c(args, '-g', 'famConstraint.tre')
	}
	
	# add thorough arguments
	args <- c(args, '-nstop', '500', '-pers', '0.2')

	if (length(tax) < 500) {
		args <- c(args, '--runs', nIter, '-pre', paste0(thoroughClades[i], '_thorough'), '-t', bestDefaultTree)
	} else {
		
		args <- c(args, '--seed')
		args <- lapply(1:nIter, function(x) c(args, paste(rep(x, 3), collapse = ''), '-pre', paste0(thoroughClades[i], '_thorough_', x), '-t', bestDefaultTree))
		
	}
	
	# runtime
	# STANDARD IQTREE - for a single run
	runTime <- "2:00:00:00" # nIter sequential runs
	
	if (length(tax) > 100) { # nIter sequential runs
		runTime <- "5:00:00:00"
	}

	if (length(tax) >= 500) { # nIter independent runs
		runTime <- "3:00:00:00"
	}

	if (length(tax) >= 1000) { # nIter independent runs
		runTime <- "4:00:00:00"
	}

	# modify PBS template
	PBScommand <- args
	
	if (!is.list(args)) {
		PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
		tmpPBS <- pbsTemplate
		tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', thoroughClades[i], grep('jobName', tmpPBS, value=TRUE))
		tmpPBS[grep('#PBS -l nodes', tmpPBS)] <- gsub('(ppn=)(\\d\\d?)', paste0('ppn=', nCores1), grep('#PBS -l nodes', tmpPBS, value = TRUE))
		tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
		tmpPBS[grep('vmem', tmpPBS)] <- gsub('(vmem=)(\\d\\d?)', paste0('vmem=', mem1), grep('vmem', tmpPBS, value=TRUE))
		tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', thoroughClades[i], ';')
		tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')

		write(tmpPBS, file = paste0(pbsDir, '/', thoroughClades[i], '_con.pbs'))
	} else {
		
		for (j in 1:nIter) {
			PBScommand[[j]][1] <- '/N/u/ptitle/Karst/iqtree2'
			
			tmpPBS <- pbsTemplate
			tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(thoroughClades[i], j), grep('jobName', tmpPBS, value=TRUE))
			tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
			tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', thoroughClades[i], ';')
			tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand[[j]], collapse = ' '), ';')

			write(tmpPBS, file = paste0(pbsDir, '/', thoroughClades[i], '_', j, '_con.pbs'))
			
		}
	}
			
	# Write PBS file for unconstrained analysis
	if (hasConstraint) {
		if (!is.list(args)) {
			PBScommand <- args[!grepl('^-g$|famConstraint.tre$', args)]
			PBScommand[grep('-pre', PBScommand) + 1] <- paste0(PBScommand[grep('-pre', PBScommand) + 1], '_uncon')
			PBScommand[grep('treefile$', PBScommand)] <- bestDefaultTree_uncon
	
			PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
			tmpPBS <- pbsTemplate
			tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(thoroughClades[i], '_uncon'), grep('jobName', tmpPBS, value=TRUE))
			tmpPBS[grep('#PBS -l nodes', tmpPBS)] <- gsub('(ppn=)(\\d\\d?)', paste0('ppn=', nCores1), grep('#PBS -l nodes', tmpPBS, value = TRUE))
			tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
			tmpPBS[grep('vmem', tmpPBS)] <- gsub('(vmem=)(\\d\\d?)', paste0('vmem=', mem1), grep('vmem', tmpPBS, value=TRUE))
			tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', thoroughClades[i], ';')
			tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')
	
			write(tmpPBS, file = paste0(pbsDir, '/', thoroughClades[i], '_uncon.pbs'))
		} else {
			
			for (j in 1:nIter) {
				PBScommand[[j]] <- args[[j]][!grepl('^-g$|famConstraint.tre$', args[[j]])]
				PBScommand[[j]][grep('-pre', PBScommand[[j]]) + 1] <- paste0(PBScommand[[j]][grep('-pre', PBScommand[[j]]) + 1], '_uncon')
				PBScommand[[j]][grep('treefile$', PBScommand[[j]])] <- bestDefaultTree_uncon 
	
				PBScommand[[j]][1] <- '/N/u/ptitle/Karst/iqtree2'
				
				tmpPBS <- pbsTemplate
				tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', paste0(thoroughClades[i], '_uncon', j), grep('jobName', tmpPBS, value=TRUE))
				tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
				tmpPBS[grep('familyPath', tmpPBS)] <- paste0('cd ./families/', thoroughClades[i], ';')
				tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand[[j]], collapse = ' '), ';')
	
				write(tmpPBS, file = paste0(pbsDir, '/', thoroughClades[i], '_', j, '_uncon.pbs'))
				
			}
		}
	}	
}

allPBS <- list.files(pbsDir, pattern = '\\.pbs$')
write(paste0('qsub ', allPBS, ';'), file = qsubFile)

##################################################################
# # Move some new files into main folders

loc1 <- '~/Dropbox/squamatePhylo/2019/topologyScripts/thoroughRuns/'
loc2 <- paste0(basedir, 'tier1/families/')

for (i in 1:length(families)) {
	
	loc1a <- paste0(loc1, families[i])
	loc2a <- paste0(loc2, families[i])
	
	f <- list.files(loc1a, full.names = TRUE)
	f <- f[grep('thorough', basename(f))]
	
	if (!any(grepl('thorough', list.files(loc2a)))) {
		file.copy(f, loc2a)
	}
}



##################################################################
# # Read in default + allnni results, generate plots, export best trees.

require(phangorn)

bestTreeDir <- paste0(basedir, 'tier1/fulltree/bestFamilyTrees')
allFamilyTreesDir <- paste0(basedir, 'tier1/allFamilyTrees')
pdfFile <- paste0(basedir, 'tier1/tier1_families_iqtree.pdf')

# thoroughPath <- '~/Dropbox/squamatePhylo/2019/topologyScripts/thoroughRuns/'
thoroughPath <- paste0(basedir, 'tier1/families/')

writeFiles <- FALSE
plot <- FALSE
writeAllTrees <- FALSE

# collect all trees and reported likelihood scores

scoresList <- vector('list', length(families))
treesList <- vector('list', length(families))
names(scoresList) <- families
names(treesList) <- families
constraintSizes <- numeric(length(families))

for (i in 1:length(families)) {
	
	tax <- unique(taxonTable[which(taxonTable$family == families[i]), 'repdbTaxon'])
	
	# are any taxa not in the concatenated alignment?
	tax <- intersect(tax, alnTaxa)
	
	# remove excluded taxa
	tax <- setdiff(tax, excludeTaxa)
		
	# subset constraint tree
	if (length(intersect(genomTree$tip.label, tax)) > 3) {
		hasConstraint <- TRUE
	} else {
		hasConstraint <- FALSE
	}
		
	message('\t', families[i], ' -- has constraint: ', hasConstraint)
	
	familyDir <- paste0(basedir, 'tier1/families/', families[i])
	
	setwd(familyDir)
	
	if (hasConstraint) {
		constraintSizes[i] <- Ntip(read.tree('famConstraint.tre'))
	} else {
		constraintSizes[i] <- 0
	}
	
	# CONSTRAINED
	# default runs
	clusterFiles <- list.files()
	clusterFiles <- clusterFiles[grep('default', clusterFiles)]
	logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
	logFiles <- logFiles[!grepl('uncon', logFiles)]
	treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
	treeFiles1 <- treeFiles1[!grepl('uncon', treeFiles1)]
	iterTreeFiles <- grep('runtrees$', clusterFiles, value = TRUE)
	iterTreeFiles <- iterTreeFiles[!grepl('uncon', iterTreeFiles)]
	
	# parallel runs
	if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
				
		# enforce ordering
		logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
		treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]
		
		defaultRuns <- numeric(50)
		for (j in 1:50) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			defaultRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		iterTreeFiles <- treeFiles1
		
	# sequential runs	
	} else {
		if (length(logFiles) > 1) stop()
		labs1 <- paste0('default_constrained_', families[i])
		tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
		zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
		zz <- do.call(rbind, strsplit(zz, '\t'))
		defaultRuns <- as.numeric(zz[,2])
	}	
	
	if (hasConstraint) {
		# UNCONSTRAINED
		# default runs
		clusterFiles <- list.files()
		clusterFiles <- clusterFiles[grep('default', clusterFiles)]
		logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
		logFiles <- logFiles[grepl('uncon', logFiles)]
		treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
		treeFiles1 <- treeFiles1[grepl('uncon', treeFiles1)]
		iterTreeFiles_uncon <- grep('runtrees$', clusterFiles, value = TRUE)
		iterTreeFiles_uncon <- iterTreeFiles_uncon[grepl('uncon', iterTreeFiles_uncon)]
		
		# parallel runs
		if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {

		# enforce ordering
		logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
		treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]

		defaultRuns_uncon <- numeric(50)
		for (j in 1:50) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			defaultRuns_uncon[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		iterTreeFiles_uncon <- treeFiles1			
			
		# sequential runs	
		} else {
			if (length(logFiles) > 1) stop()
			labs1 <- paste0('default_constrained_', families[i])
			tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
			zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
			zz <- do.call(rbind, strsplit(zz, '\t'))
			defaultRuns_uncon <- as.numeric(zz[,2])	
		}
	} else {
		defaultRuns_uncon <- NA
	}
		
	# ALLNNI-ADD-ON runs
	
	# CONSTRAINED
	clusterFiles <- list.files()
	clusterFiles <- clusterFiles[grep('allnni', clusterFiles)]
	logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
	logFiles <- logFiles[!grepl('uncon', logFiles)]
	logFiles <- logFiles[!grepl('allnn?iOnly', logFiles)]
	treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
	treeFiles1 <- treeFiles1[!grepl('uncon', treeFiles1)]
	treeFiles1 <- treeFiles1[!grepl('allnn?iOnly', treeFiles1)]
	iterTreeFiles_allnni <- grep('runtrees$', clusterFiles, value = TRUE)
	iterTreeFiles_allnni <- iterTreeFiles_allnni[!grepl('uncon|allnn?iOnly', iterTreeFiles_allnni)]

	# parallel runs
	if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {

		# enforce ordering
		logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
		treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]
		
		allnniRuns <- numeric(10)
		for (j in 1:10) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			allnniRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		iterTreeFiles_allnni <- treeFiles1
		
	# sequential runs	
	} else {
		if (length(logFiles) > 1) stop()
		labs1 <- paste0('default_constrained_', families[i])
		tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
		zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
		zz <- do.call(rbind, strsplit(zz, '\t'))
		allnniRuns <- as.numeric(zz[,2])
	}	
	
	if (hasConstraint) {
		# UNCONSTRAINED
		# allnni runs
		clusterFiles <- list.files()
		clusterFiles <- clusterFiles[grep('allnni', clusterFiles)]
		logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
		logFiles <- logFiles[grepl('uncon', logFiles)]
		logFiles <- logFiles[!grepl('allnn?iOnly', logFiles)]
		treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
		treeFiles1 <- treeFiles1[grepl('uncon', treeFiles1)]
		treeFiles1 <- treeFiles1[!grepl('allnn?iOnly', treeFiles1)]
		iterTreeFiles_allnni_uncon <- grep('runtrees$', clusterFiles, value = TRUE)
		iterTreeFiles_allnni_uncon <- iterTreeFiles_allnni_uncon[grepl('uncon|allnn?iOnly', iterTreeFiles_allnni_uncon)]
		
		# parallel runs
		if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {

		# enforce ordering
		logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
		treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]

		allnniRuns_uncon <- numeric(10)
		for (j in 1:10) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			allnniRuns_uncon[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		iterTreeFiles_allnni_uncon <- treeFiles1
			
		# sequential runs	
		} else {
			if (length(logFiles) > 1) stop()
			labs1 <- paste0('default_constrained_', families[i])
			tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
			zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
			zz <- do.call(rbind, strsplit(zz, '\t'))
			allnniRuns_uncon <- as.numeric(zz[,2])	
		}
	} else {
		allnniRuns_uncon <- NA
	}
	
	# ADDITIONAL THOROUGH RUNS
	famThoroughPath <- paste0(thoroughPath, families[i])
	# if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
		# famThoroughPath <- getwd()
	# }
	if (any(grepl('thorough', list.files(famThoroughPath)))) {
	
		# constrained
		clusterFiles <- list.files(famThoroughPath, full.names = TRUE)
		clusterFiles <- clusterFiles[grep('thorough', basename(clusterFiles))]
		logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
		logFiles <- logFiles[!grepl('uncon', logFiles)]
		treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
		treeFiles1 <- treeFiles1[!grepl('uncon', treeFiles1)]
		iterTreeFiles_thorough <- grep('runtrees$', clusterFiles, value = TRUE)
		iterTreeFiles_thorough <- iterTreeFiles_thorough[!grepl('uncon', iterTreeFiles_thorough)]

		# parallel runs
		if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
	
			# enforce ordering
			logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
			treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]
			
			thoroughRuns <- numeric(length(logFiles))
			for (j in 1:length(logFiles)) {
				tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
				zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
				thoroughRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
			}
			iterTreeFiles_thorough <- treeFiles1
	
		} else {
			if (length(logFiles) > 1) stop()
			labs1 <- paste0('thorough_constrained_', families[i])
			tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
			zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
			zz <- do.call(rbind, strsplit(zz, '\t'))
			thoroughRuns <- as.numeric(zz[,2])	
		}

		if (hasConstraint) {
			# unconstrained
			clusterFiles <- list.files(famThoroughPath, full.names = TRUE)
			clusterFiles <- clusterFiles[grep('thorough', basename(clusterFiles))]
			logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
			logFiles <- logFiles[grepl('uncon', logFiles)]
			treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
			treeFiles1 <- treeFiles1[grepl('uncon', treeFiles1)]
			iterTreeFiles_thorough_uncon <- grep('runtrees$', clusterFiles, value = TRUE)
			iterTreeFiles_thorough_uncon <- iterTreeFiles_thorough_uncon[grepl('uncon', iterTreeFiles_thorough_uncon)]
	
			# parallel runs
			if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
		
				# enforce ordering
				logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
				treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]
				
				thoroughRuns_uncon <- numeric(length(logFiles))
				for (j in 1:length(logFiles)) {
					tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
					zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
					thoroughRuns_uncon[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
				}
				iterTreeFiles_thorough_uncon <- treeFiles1
		
			} else {
				if (length(logFiles) > 1) stop()
				labs1 <- paste0('thorough_unconstrained_', families[i])
				tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
				zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
				zz <- do.call(rbind, strsplit(zz, '\t'))
				thoroughRuns_uncon <- as.numeric(zz[,2])	
			}
		} else {
			thoroughRuns_uncon <- NA
		}
	} else {
		thoroughRuns <- NA
		thoroughRuns_uncon <- NA
	}
	
	# ALLNNI-ONLY runs
	
	# CONSTRAINED
	clusterFiles <- list.files()
	clusterFiles <- clusterFiles[grep('allnn?iOnly', clusterFiles)]
	logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
	logFiles <- logFiles[!grepl('uncon', logFiles)]
	treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
	treeFiles1 <- treeFiles1[!grepl('uncon', treeFiles1)]
	iterTreeFiles_allnniOnly <- grep('runtrees$', clusterFiles, value = TRUE)
	iterTreeFiles_allnniOnly <- iterTreeFiles_allnniOnly[!grepl('uncon', iterTreeFiles_allnniOnly)]

	# parallel runs
	if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {

		# enforce ordering
		logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
		treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]
		
		allnniOnlyRuns <- numeric(5)
		for (j in 1:5) {
			tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
			zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
			allnniOnlyRuns[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
		}
		iterTreeFiles_allnniOnly <- treeFiles1
		
	# sequential runs	
	} else {
		if (length(logFiles) > 1) stop()
		labs1 <- paste0('default_constrained_', families[i])
		tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
		zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
		zz <- do.call(rbind, strsplit(zz, '\t'))
		allnniOnlyRuns <- as.numeric(zz[,2])
	}	
	
	if (hasConstraint) {
		# UNCONSTRAINED
		# allnni runs
		clusterFiles <- list.files()
		clusterFiles <- clusterFiles[grep('allnn?iOnly', clusterFiles)]
		logFiles <- grep('iqtree$', clusterFiles, value = TRUE)
		logFiles <- logFiles[grepl('uncon', logFiles)]
		treeFiles1 <- grep('treefile$', clusterFiles, value = TRUE)
		treeFiles1 <- treeFiles1[grepl('uncon', treeFiles1)]
		iterTreeFiles_allnniOnly_uncon <- grep('runtrees$', clusterFiles, value = TRUE)
		iterTreeFiles_allnniOnly_uncon <- iterTreeFiles_allnniOnly_uncon[grepl('uncon', iterTreeFiles_allnniOnly_uncon)]
		
		# parallel runs
		if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {

			# enforce ordering
			logFiles <- logFiles[order(as.numeric(stringr::str_extract(logFiles, '\\d\\d?')))]
			treeFiles1 <- treeFiles1[order(as.numeric(stringr::str_extract(treeFiles1, '\\d\\d?')))]
	
			allnniOnlyRuns_uncon <- numeric(5)
			for (j in 1:5) {
				tmp <- scan(logFiles[j], what = 'character', sep = '\n', quiet = TRUE)
				zz <- grep('Log-likelihood of the tree:', tmp, value = TRUE)
				allnniOnlyRuns_uncon[j] <- as.numeric(gsub("(Log-likelihood of the tree: )(-\\d+\\.\\d+)(\\s\\(.+)", "\\2", zz))
			}
			iterTreeFiles_allnniOnly_uncon <- treeFiles1
			
		# sequential runs	
		} else {
			if (length(logFiles) > 1) stop()
			labs1 <- paste0('default_constrained_', families[i])
			tmp <- scan(logFiles, what = 'character', sep = '\n', quiet = TRUE)
			zz <- tmp[(grep('MULTIPLE RUNS', tmp) + 3):(grep('MAXIMUM LIKELIHOOD TREE', tmp) - 1)]
			zz <- do.call(rbind, strsplit(zz, '\t'))
			allnniOnlyRuns_uncon <- as.numeric(zz[,2])	
		}
	} else {
		allnniOnlyRuns_uncon <- NA
	}
	
	##############	
	# pull together all scores and trees
	
	# constrained
	loglik_con <- c(defaultRuns, allnniRuns, thoroughRuns, allnniOnlyRuns)
	#loglik_con <- loglik_con[!is.na(loglik_con)]
	if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
		alltrees_con <- c(as.list(iterTreeFiles), as.list(iterTreeFiles_allnni), as.list(iterTreeFiles_thorough), as.list(iterTreeFiles_allnniOnly))
		alltrees_con <- lapply(alltrees_con, read.tree)
	} else {
		alltrees_con <- c(read.tree(iterTreeFiles), read.tree(iterTreeFiles_allnni), read.tree(iterTreeFiles_thorough), read.tree(iterTreeFiles_allnniOnly))
	}
	
	# unconstrained
	if (hasConstraint) {
		loglik_uncon <- c(defaultRuns_uncon, allnniRuns_uncon, thoroughRuns_uncon, allnniOnlyRuns_uncon)
		#loglik_uncon <- loglik_uncon[!is.na(loglik_uncon)]
		if (families[i] %in% c('Colubridae', 'Scincidae', 'Gekkonidae')) {
			alltrees_uncon <- c(as.list(iterTreeFiles_uncon), as.list(iterTreeFiles_allnni_uncon), as.list(iterTreeFiles_thorough_uncon), as.list(iterTreeFiles_allnniOnly_uncon))
			alltrees_uncon <- lapply(alltrees_uncon, read.tree)
		} else {
			alltrees_uncon <- c(read.tree(iterTreeFiles_uncon), read.tree(iterTreeFiles_allnni_uncon), read.tree(iterTreeFiles_thorough_uncon), read.tree(iterTreeFiles_allnniOnly_uncon))
		}
	}
	
	# store all scores and trees

	if (hasConstraint) {
		famScores_con <- data.frame(family = families[i], constrained = TRUE, type = NA, score = loglik_con)
		famScores_con$type = c(rep('default', length(defaultRuns)), rep('allnni', length(allnniRuns)), rep('thorough', length(thoroughRuns)), rep('allnniOnly', length(allnniOnlyRuns)))

		famScores_uncon <- data.frame(family = families[i], constrained = FALSE, type = NA, score = loglik_uncon)
		famScores_uncon$type = c(rep('default', length(defaultRuns_uncon)), rep('allnni', length(allnniRuns_uncon)), rep('thorough', length(thoroughRuns_uncon)), rep('allnniOnly', length(allnniOnlyRuns_uncon)))
		
		max(sapply(split(famScores_con$score, famScores_con$type), max)) - sapply(split(famScores_con$score, famScores_con$type), max)
		max(sapply(split(famScores_uncon$score, famScores_uncon$type), max)) - sapply(split(famScores_uncon$score, famScores_uncon$type), max)
		
		treesList[[i]] <- list(con = alltrees_con, uncon = alltrees_uncon)
		
		
	} else {
		famScores_con <- data.frame(family = families[i], constrained = TRUE, type = NA, score = rep(NA, length(loglik_con)))
		famScores_con$type = c(rep('default', length(defaultRuns)), rep('allnni', length(allnniRuns)), rep('thorough', length(thoroughRuns)), rep('allnniOnly', length(allnniOnlyRuns)))

		famScores_uncon <- data.frame(family = families[i], constrained = FALSE, type = NA, score = loglik_con)
		famScores_uncon$type = c(rep('default', length(defaultRuns)), rep('allnni', length(allnniRuns)), rep('thorough', length(thoroughRuns)), rep('allnniOnly', length(allnniOnlyRuns)))
		
		treesList[[i]] <- list(con = NA, uncon = alltrees_con)
		
	}
	
	scoresList[[i]] <- list(famScores_con, famScores_uncon)
	
	# clear these variables, to avoid any unintended confusion
	rm(loglik_con, defaultRuns, allnniRuns, thoroughRuns, allnniOnlyRuns, loglik_uncon, defaultRuns_uncon, allnniRuns_uncon, thoroughRuns_uncon, allnniOnlyRuns_uncon, alltrees_con, iterTreeFiles, iterTreeFiles_allnni, iterTreeFiles_thorough, iterTreeFiles_allnniOnly, alltrees_uncon, iterTreeFiles_uncon, iterTreeFiles_allnni_uncon, iterTreeFiles_thorough_uncon, iterTreeFiles_allnniOnly_uncon)	
}
	
	
sapply(scoresList, function(x) length(na.omit(do.call(rbind, x)[,4])))
sapply(treesList, function(x) length(unlist(x, recursive = FALSE)))

###########################
# write trees and scores
for (i in 1:length(families)) {
	
	# con
	tmp <- treesList[[i]]$con
	if (!anyNA(tmp)) {
		class(tmp) <- 'multiPhylo'
		write.tree(tmp, file = paste0(allFamilyTreesDir, '/', families[i], '_con.trees'))
	}
	
	# uncon
	tmp <- treesList[[i]]$uncon
	class(tmp) <- 'multiPhylo'
	write.tree(tmp, file = paste0(allFamilyTreesDir, '/', families[i], '_uncon.trees'))
}
	
saveRDS(scoresList, paste0(allFamilyTreesDir, '/reportedLikelihoods.rds'))
saveRDS(constraintSizes, paste0(allFamilyTreesDir, 'constraintSizes.rds'))
	
#############################
# Generate plots

scoresList <- readRDS(paste0(allFamilyTreesDir, '/reportedLikelihoods.rds'))
constraintSizes <- readRDS(paste0(allFamilyTreesDir, 'constraintSizes.rds'))

families <- names(scoresList)

treesList <- list()
for (i in 1:length(families)) {
	f <- list.files(allFamilyTreesDir, pattern = families[i], full.name = TRUE)
	f_uncon <- f[grep('uncon', basename(f))]
	f_con <- f[grep('_con', basename(f))]
	f_con <- ifelse(length(f_con) == 0, NA, f_con)
	if (!is.na(f_con)) {
		f_con <- read.tree(f_con)
	}
	f_uncon <- read.tree(f_uncon)
	treesList[[i]] <- list(con = f_con, uncon = f_uncon)
}

pdf(pdfFile, width = 12, height = 7, onefile = TRUE)
	
for (i in 1:length(families)) {
	
	message(i)

	nConTips <- constraintSizes[i]
	
	tmp <- scoresList[[i]] # con first, uncon second
	loglik_con <- tmp[[1]]$score
	loglik_uncon <- tmp[[2]]$score
	
	contrees <- treesList[[i]]$con
	uncontrees <- treesList[[i]]$uncon
	
	if (all(is.na(loglik_con))) {
		hasConstraint <- FALSE
		rfdists_con <- NA
		rfdists_uncon <- sapply(uncontrees, function(x) RF.dist(uncontrees[[which.max(loglik_uncon)]], x, normalize = FALSE))
	} else {
		hasConstraint <- TRUE
		rfdists_con <- sapply(contrees, function(x) RF.dist(contrees[[which.max(loglik_con)]], x, normalize = FALSE))
		rfdists_uncon <- sapply(uncontrees, function(x) RF.dist(uncontrees[[which.max(loglik_uncon)]], x, normalize = FALSE))
	}
	
	par(mfrow = c(1,2))
	
	plot.new()
	plot.window(xlim = c(1, length(loglik_con)), ylim = range(c(loglik_con, loglik_uncon), na.rm = TRUE))
	axlabels <- 1:max(length(loglik_con), length(loglik_uncon))
	#axlabels <- c(1:50, 1:10, 1:length(thoroughRuns))[1:length(loglik_con)]
	#axis(1, at = 1:length(loglik_con), labels = axlabels, las=2, cex.axis=0.7)
	axis(1, at = axlabels, labels = axlabels, las=2, cex.axis=0.7)
	axis(2, cex.axis = 0.7)
	box()

	points(1:length(loglik_uncon), loglik_uncon, col = gray(0.5), pch = 2)
	points(which.max(loglik_uncon), loglik_uncon[which.max(loglik_uncon)], pch = 17, col = gray(0.5))
	abline(h = loglik_uncon[which.max(loglik_uncon)], lwd = 0.3, lty=3, col = gray(0.5))

	if (hasConstraint) {
		points(1:length(loglik_con), loglik_con)
		points(which.max(loglik_con), loglik_con[which.max(loglik_con)], pch = 20)
		abline(h = loglik_con[which.max(loglik_con)], lwd = 0.3, lty=3)
	}

	abline(v = c(50.5, 60.5, ifelse(length(loglik_con) == 70, 65.5, 66.5)), lty=3)		

	mtext('likelihood', side = 2, line = 2.2, cex = 0.75)
	legend(x = grconvertX(0, from = 'npc'), y = grconvertY(1.1, from = 'npc'), legend = c('con', 'uncon'), col = c('black', 'gray50'), pch = c(1, 2), bty = 'n', xpd = NA)
	mtext(paste0('best uncon lk: ', round(max(loglik_uncon), 2)), side = 1, line = 2, cex=0.7)

	mtext(paste0('best con lk: ', round(max(loglik_con), 2)), side = 1, line = 2.7, cex=0.7)
	title(main = 'likelihood')
	
	plot.new()
	plot.window(xlim = c(1, length(rfdists_uncon)), ylim = rev(range(c(rfdists_con , rfdists_uncon), na.rm = TRUE)))
	#axlabels <- c(1:50, 1:10, 1:length(thoroughRuns))[1:length(loglik_con)]
	axlabels <- 1:max(length(loglik_uncon), length(loglik_uncon))
	axis(1, at = 1:length(rfdists_uncon), labels = axlabels, las=2, cex.axis=0.7)
	axis(2, cex.axis = 0.7)
	box()

	points(1:length(rfdists_uncon), rfdists_uncon, col = gray(0.5), pch = 2)
	points(which.max(loglik_uncon), rfdists_uncon[which.max(loglik_uncon)], , pch = 17, col = gray(0.5))

	if (hasConstraint) {
		points(1:length(rfdists_con), rfdists_con)
		points(which.max(loglik_con), rfdists_con[which.max(loglik_con)], pch = 20)
	}

	# abline(v = 10.5, lty=3)	
	mtext('RF distance', side = 2, line = 2.2, cex = 0.75)
	# mtext(families[i], outer = TRUE, adj = 0.5, line = -3, font = 2)
	
	abline(v = c(50.5, 60.5, ifelse(length(loglik_con) == 70, 65.5, 66.5)), lty=3)
	title(main = 'RF dist')
	
	mtext(paste0(families[i], ' (', Ntip(uncontrees[[1]]) - length(outgroupTaxa), ')'), outer = TRUE, line = -3, font = 2, cex = 1.2)
	mtext(paste0('n constrained tips: ', nConTips), side = 1, outer = TRUE, line = -3, font = 2, cex = 1)
}

dev.off()
	

	
# Prep script to evaluate likelihoods using raxml-ng

allcmd <- c()
for (i in 1:length(families)) {
	
	f <- list.files(allFamilyTreesDir, pattern = families[i])
	f_uncon <- f[grep('uncon', f)]
	f_con <- f[grep('_con', f)]
	f_con <- ifelse(length(f_con) == 0, NA, f_con)
	
	# ~/raxml-ng --evaluate --msa famConcat.aln --model ../concatenated.partitions.raxml --threads 12 --prefix evalRaxmlng_family --opt-model on --opt-branches off --tree trees.trees
	
	if (!is.na(f_con)) {
		
		cmd_con <- paste('~/raxml-ng', '--evaluate', '--msa', paste0('../', families[i], '/famConcat.aln'), '--model', '../concatenated.partitions.raxml', '--threads', 16, '--opt-model', 'on', '--opt-branches', 'off', '--prefix', paste0('evalRaxmlng_', families[i], '_con'), '--tree', paste0('../famTrees/', f_con))
		
		allcmd <- c(allcmd, cmd_con)
		
	}
	
	cmd_uncon <- paste('~/raxml-ng', '--evaluate', '--msa', paste0('../', families[i], '/famConcat.aln'), '--model', '../concatenated.partitions.raxml', '--threads', 16, '--opt-model', 'on', '--opt-branches', 'off', '--prefix', paste0('evalRaxmlng_', families[i], '_uncon'), '--tree', paste0('../famTrees/', f_uncon))
	
	allcmd <- c(allcmd, cmd_uncon)
}

allcmd <- paste0(allcmd, ';')

slurmHeader <- 
'#!/bin/bash
#
#SBATCH --job-name=evalAllTrees
#SBATCH --output=evalAllTrees.txt
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=167:00:00
#SBATCH -p extended-28core
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=pascal.title@stonybrook.edu

cd /gpfs/scratch/ptitle/eval/raxmlEval;

'

slurmFile <- paste0(basedir, 'tier1/allFamilyTrees/raxml_Eval.slurm')

write(slurmHeader, file = slurmFile)
write(allcmd, file = slurmFile, append = TRUE)

#################################
# Read in raxml-ng likelihoods

setwd('~/Dropbox/squamatePhylo/2019/topologyScripts/raxmlEval')

allFamilyTreesDir <- paste0(basedir, 'tier1/allFamilyTrees')

scoresList <- readRDS(paste0(allFamilyTreesDir, '/reportedLikelihoods.rds'))

for (i in 1:length(families)) {
	
	message('\t', families[i])
	
	f <- list.files(pattern = families[i])
	f <- grep('\\.log$', f, value = TRUE)
	f_uncon <- f[grep('uncon', f)]
	f_con <- f[grep('_con', f)]
	f_con <- ifelse(length(f_con) == 0, NA, f_con)
	
	f_uncon <- scan(f_uncon, what = 'character', sep = '\n', quiet = TRUE)
	f_uncon <- grep('\\#\\d\\d?, final logLikelihood', f_uncon, value = TRUE)
	f_uncon <- stringr::str_extract(f_uncon, '-\\d+\\.\\d+')
	f_uncon <- as.numeric(f_uncon)
	#scoresList[[i]][[2]]$raxml <- f_uncon
	scoresList[[i]][[2]] <- cbind.data.frame(scoresList[[i]][[2]], raxml = f_uncon)
	
	if (!is.na(f_con)) {
		
		f_con <- scan(f_con, what = 'character', sep = '\n', quiet = TRUE)
		f_con <- grep('\\#\\d\\d?, final logLikelihood', f_con, value = TRUE)
		f_con <- stringr::str_extract(f_con, '-\\d+\\.\\d+')
		f_con <- as.numeric(f_con)
		#scoresList[[i]][[1]]$raxml <- f_con
		scoresList[[i]][[1]] <- cbind.data.frame(scoresList[[i]][[1]], raxml = f_con)		
	}
}


# count monophyletic genera
require(MonoPhy)

tax <- unique(taxonTable[, c('repdbTaxon', 'genus', 'family')])
colnames(tax) <- c('species', 'genus', 'family')
rownames(tax) <- tax$species

treesList <- list()
for (i in 1:length(families)) {
	f <- list.files(allFamilyTreesDir, pattern = families[i], full.name = TRUE)
	f_uncon <- f[grep('uncon', basename(f))]
	f_con <- f[grep('_con', basename(f))]
	f_con <- ifelse(length(f_con) == 0, NA, f_con)
	if (!is.na(f_con)) {
		f_con <- read.tree(f_con)
	}
	f_uncon <- read.tree(f_uncon)
	treesList[[i]] <- list(con = f_con, uncon = f_uncon)
}


for (i in 1:length(families)) {
	
	message('\t', families[i])
	
	famTaxa <- tax[which(tax$family == families[i]),]

	unconMonoGenera <- integer(length(treesList[[i]]$uncon))
	for (j in 1:length(treesList[[i]]$uncon)) {
		tmp <- root(treesList[[i]]$uncon[[j]], 'Sphenodon_punctatus', resolve.root = TRUE)
		tmp <- keep.tip(tmp, intersect(tmp$tip.label, famTaxa$species))
		unconMonoGenera[j] <- as.numeric(GetSummaryMonophyly(AssessMonophyly(tmp, famTaxa[tmp$tip.label,]))$genus['Monophyletic', 'Taxa'])
	}
	unconMonoGenera[is.na(unconMonoGenera)] <- 0
	
	scoresList[[i]][[2]] <- cbind.data.frame(scoresList[[i]][[2]], monoGenera = unconMonoGenera)
	
	if (!anyNA(treesList[[i]]$con)) {
	
		conMonoGenera <- integer(length(treesList[[i]]$con))
		for (j in 1:length(treesList[[i]]$con)) {
			tmp <- root(treesList[[i]]$con[[j]], 'Sphenodon_punctatus', resolve.root = TRUE)
			tmp <- keep.tip(tmp, intersect(tmp$tip.label, famTaxa$species))
			conMonoGenera[j] <- as.numeric(GetSummaryMonophyly(AssessMonophyly(tmp, famTaxa[tmp$tip.label,]))$genus['Monophyletic', 'Taxa'])
		}
		conMonoGenera[is.na(conMonoGenera)] <- 0
		scoresList[[i]][[1]] <- cbind.data.frame(scoresList[[i]][[1]], monoGenera = conMonoGenera)	
	}
}

for (i in 1:length(families)) {
	
	if (!anyNA(scoresList[[i]][[1]]$score)) {
		if (anyNA(scoresList[[i]][[1]][,6])) {
			scoresList[[i]][[1]][,6][is.na(scoresList[[i]][[1]][,6])] <- 0
		}
		
		
	}
	
	if (anyNA(scoresList[[i]][[2]][,6])) {
		scoresList[[i]][[2]][,6][is.na(scoresList[[i]][[2]][,6])] <- 0
	}
	
}

pdf(paste0(basedir, 'tier1/tier1_families_loglikCompared.pdf'), width = 10, height = 10, onefile = TRUE)

for (i in 1:length(families)) {

	par(mfrow=c(2,2))
	
	plot.new()
	if (!is.na(scoresList[[i]][[1]][1, 'score'])) {
		likRange <- range(c(scoresList[[i]][[1]]$score, scoresList[[i]][[1]]$raxml))
		plot.window(xlim = likRange, ylim = likRange)
		axis(1, cex.axis = 0.5)
		axis(2, cex.axis = 0.5)
		box()
		points(scoresList[[i]][[1]]$score, scoresList[[i]][[1]]$raxml)
		abline(a = 0, b = 1, lty = 3)
		mtext('IQTREE reported scores', side = 1, line = 2.5)
		mtext('RAXML-NG scores', side = 2, line = 2.5)
		mtext('constrained', side = 3, line = 1)
	} 	
	
	plot.new()
	likRange <- range(c(scoresList[[i]][[2]]$score, scoresList[[i]][[2]]$raxml))
	plot.window(xlim = likRange, ylim = likRange)
	axis(1, cex.axis = 0.5)
	axis(2, cex.axis = 0.5)
	box()
	points(scoresList[[i]][[2]]$score, scoresList[[i]][[2]]$raxml)
	abline(a = 0, b = 1, lty = 3)
	mtext('IQTREE reported scores', side = 1, line = 2.5)
	mtext('RAXML-NG scores', side = 2, line = 2.5)
	mtext('unconstrained', side = 3, line = 1)
	mtext(families[i], side = 3, outer = TRUE, line = -2, font = 2)
	
	# constrained zoom-in
	plot.new()
	if (!is.na(scoresList[[i]][[1]][1, 'score'])) {
		plot.window(xlim = range(scoresList[[i]][[1]]$score), ylim = range(scoresList[[i]][[1]]$raxml))
		axis(1, cex.axis = 0.5)
		axis(2, cex.axis = 0.5)
		box()
		points(scoresList[[i]][[1]]$score, scoresList[[i]][[1]]$raxml)
		legend('bottomright', legend = 'zoom-in', border = NULL, bty = 'n')
		#abline(lm(scoresList[[i]][[1]]$raxml ~ scoresList[[i]][[1]]$score), lty = 3)
		mtext('IQTREE reported scores', side = 1, line = 2.5)
		mtext('RAXML-NG scores', side = 2, line = 2.5)
		mtext('constrained', side = 3, line = 1)
	} 	
	
	# unconstrained zoom-in
	plot.new()
	plot.window(xlim = range(scoresList[[i]][[2]]$score), ylim = range(scoresList[[i]][[2]]$raxml))
	axis(1, cex.axis = 0.5)
	axis(2, cex.axis = 0.5)
	box()
	points(scoresList[[i]][[2]]$score, scoresList[[i]][[2]]$raxml)
	legend('bottomright', legend = 'zoom-in', border = NULL, bty = 'n')
	#abline(lm(scoresList[[i]][[2]]$raxml ~ scoresList[[i]][[2]]$score), lty = 3)
	mtext('IQTREE reported scores', side = 1, line = 2.5)
	mtext('RAXML-NG scores', side = 2, line = 2.5)
	mtext('unconstrained', side = 3, line = 1)

}

dev.off()


# plot likelihood comparison, and relationship to monophyly

pdf(paste0(basedir, 'tier1/tier1_families_genusMonophyly_vs_loglik.pdf'), width = 10, height = 10, onefile = TRUE)

for (i in 1:length(families)) {
	
	if (!all(is.na(treesList[[i]]$con))) {
	con1 <- treesList[[i]]$con[[which.max(scoresList[[i]][[1]]$score)]]
	con2 <- treesList[[i]]$con[[which.max(scoresList[[i]][[1]]$raxml)]]
	RFdiff1 <- RF.dist(con1, con2)
	} else {
		RFdiff1 <- 0
	}
	uncon1 <- treesList[[i]]$uncon[[which.max(scoresList[[i]][[2]]$score)]]
	uncon2 <- treesList[[i]]$uncon[[which.max(scoresList[[i]][[2]]$raxml)]]
	RFdiff2 <- RF.dist(uncon1, uncon2)
	
	if (RFdiff1 > 0 | RFdiff2 > 0) {

	
		par(mfrow=c(2,2))
		
		if (!all(is.na(scoresList[[i]][[1]]$score))) {
			indIQ <- which.max(scoresList[[i]][[1]]$score)
			indRax <- which.max(scoresList[[i]][[1]]$raxml)
			plot(scoresList[[i]][[1]]$score, scoresList[[i]][[1]]$monoGenera, xlab = 'iqtree loglik', ylab = 'n mono genera')
			points(scoresList[[i]][[1]]$score[indIQ], scoresList[[i]][[1]]$monoGenera[indIQ], cex=2, col = 'blue', lwd=2, pch = 3)
			points(scoresList[[i]][[1]]$score[indRax], scoresList[[i]][[1]]$monoGenera[indRax], cex=2, col = 'orange', lwd=2, pch = 3)
			title(main = 'loglik as reported in iqtree log')
		
			plot(scoresList[[i]][[1]]$raxml, scoresList[[i]][[1]]$monoGenera, xlab = 'raxml loglik', ylab = 'n mono genera')
			points(scoresList[[i]][[1]]$raxml[indRax], scoresList[[i]][[1]]$monoGenera[indRax], cex=2, col = 'orange', lwd=2, pch = 3)
			points(scoresList[[i]][[1]]$raxml[indIQ], scoresList[[i]][[1]]$monoGenera[indIQ], cex=2, col = 'blue', lwd=2, pch = 3)
			title(main = 'loglik as reported raxml-ng')
			mtext('constrained', side = 3, outer = TRUE, at = 0.5, line = -2, font = 2)
			
			# are trees different?
			con1 <- treesList[[i]]$con[[which.max(scoresList[[i]][[1]]$score)]]
			con2 <- treesList[[i]]$con[[which.max(scoresList[[i]][[1]]$raxml)]]
			RFdiff <- RF.dist(con1, con2)
			mtext(paste0('RF dist: ', RFdiff), side = 3, at = 0.5, outer = T, line = -4)
			
		} else {
			plot.new()
			plot.new()
		}
	
		indIQ <- which.max(scoresList[[i]][[2]]$score)
		indRax <- which.max(scoresList[[i]][[2]]$raxml)
		plot(scoresList[[i]][[2]]$score, scoresList[[i]][[2]]$monoGenera, xlab = 'iqtree loglik', ylab = 'n mono genera')
		points(scoresList[[i]][[2]]$score[indIQ], scoresList[[i]][[2]]$monoGenera[indIQ], cex=2, col = 'blue', lwd=2, pch = 3)
		points(scoresList[[i]][[2]]$score[indRax], scoresList[[i]][[2]]$monoGenera[indRax], cex=2, col = 'orange', lwd=2, pch = 3)
		title(main = 'loglik as reported in iqtree log')
		legend('topleft', legend = c('iqtree', 'raxml'), col = c('blue', 'orange'), pch = 3, lwd = 2, bty = 'n')
	
		plot(scoresList[[i]][[2]]$raxml, scoresList[[i]][[2]]$monoGenera, xlab = 'raxml loglik', ylab = 'n mono genera')
		points(scoresList[[i]][[2]]$raxml[indRax], scoresList[[i]][[2]]$monoGenera[indRax], cex=2, col = 'orange', lwd=2, pch = 3)
		points(scoresList[[i]][[2]]$raxml[indIQ], scoresList[[i]][[2]]$monoGenera[indIQ], cex=2, col = 'blue', lwd=2, pch = 3)
		title(main = 'loglik as reported raxml-ng')
		mtext('unconstrained', side = 3, outer = TRUE, at = 0.5, line = -31, font = 2)
		mtext(families[i], side = 3, outer = TRUE, at = 0.1, line = -30, font = 2, cex=1.2)
		
		# are trees different?
		uncon1 <- treesList[[i]]$uncon[[which.max(scoresList[[i]][[2]]$score)]]
		uncon2 <- treesList[[i]]$uncon[[which.max(scoresList[[i]][[2]]$raxml)]]
		RFdiff <- RF.dist(uncon1, uncon2)
		mtext(paste0('RF dist: ', RFdiff), side = 3, at = 0.5, outer = T, line = -32)
	}
	
}
dev.off()


## -------------------------------------------

# Identify best tree and move to tier1/fulltree/bestFamilyTrees

# did a particular type of run produce the best trees?

# constrained
table(unlist(lapply(scoresList, function(x) x[[1]][which.max(x[[1]]$raxml), 'type'])))
table(unlist(lapply(scoresList, function(x) x[[1]][which.max(x[[1]]$score), 'type'])))

# unconstrained
table(unlist(lapply(scoresList, function(x) x[[2]][which.max(x[[2]]$raxml), 'type'])))
table(unlist(lapply(scoresList, function(x) x[[2]][which.max(x[[2]]$score), 'type'])))

for (i in 1:length(families)) {
	
	if (!all(is.na(treesList[[i]]$con))) {
		bestConTree <- treesList[[i]]$con[[which.max(scoresList[[i]][[1]]$raxml)]]
		write.tree(bestConTree, file = paste0(bestTreeDir, '/', families[i], '_con.tre'))
	} 
	
	bestUnconTree <- treesList[[i]]$uncon[[which.max(scoresList[[i]][[2]]$raxml)]]
	write.tree(bestUnconTree, file = paste0(bestTreeDir, '/', families[i], '_uncon.tre'))
}





	
	
	
##################################################################
# 
# Prepare 130-tip tree for MCMCTREE dating analyses

# In order to have a topology that will be most congruent with our genbank tree, we will create a supermatrix subset for the pre-selected 130 taxa for the time calibration tree. As with tier1 analyses, there will be a genomic topological constraint. 


# 3 locus sets for time calibration

seqFile3 <- '~/Dropbox/Oz_Crown_Ages/SortaDate/sorta_concatenated_v2/concat_n3.phy'
seqFile5 <- '~/Dropbox/Oz_Crown_Ages/SortaDate/sorta_concatenated_v2/concat_n5.phy'
seqFile10 <- '~/Dropbox/Oz_Crown_Ages/SortaDate/sorta_concatenated_v2/concat_n10.phy'

# sequence data: process and get genus_species from sample names
# seq data not needed here. Just keep sample names.

# split out sample names
seqNames <- vector('list', 3)

# 3 loci
seq <- scan(seqFile3, what = 'character', sep = '\n')
# first line is number of species/sites
seq <- seq[2:length(seq)]
for (i in 1:length(seq)) {	
	tmp <- strsplit(seq[[i]], split = '\t')[[1]]
	seqNames[[1]][i] <- tmp[[1]]
}

# 5 loci
seq <- scan(seqFile5, what = 'character', sep = '\n')
# first line is number of species/sites
seq <- seq[2:length(seq)]
for (i in 1:length(seq)) {	
	tmp <- strsplit(seq[[i]], split = '\t')[[1]]
	seqNames[[2]][i] <- tmp[[1]]
}

# 10 loci
seq <- scan(seqFile10, what = 'character', sep = '\n')
# first line is number of species/sites
seq <- seq[2:length(seq)]
for (i in 1:length(seq)) {	
	tmp <- strsplit(seq[[i]], split = '\t')[[1]]
	seqNames[[3]][i] <- tmp[[1]]
}

# merge
seqNames <- Reduce(union, seqNames)

# remove outgroups
genomOutgroups <- c('taeGut2', 'ISIS373002', 'UMFS-10956c', 'allMis1', 'chrPic1', 'galGal5', 'hg38', 'H20145a')
genomOutgroups %in% seqNames
seqNames <- setdiff(seqNames, genomOutgroups)

taxonTable <- read.csv('~/Dropbox/Oz_Crown_Ages/sqCL_taxa_repDBmatching.csv', stringsAsFactors=FALSE)
# manual matching supercedes repDB auto matching, unless no manual match found
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'repDBtaxon']
head(taxonTable)
table(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch))
taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'species']

table(seqNames %in% taxonTable$sample)
# drop the few that are not matching
seqNames <- intersect(seqNames, taxonTable$sample)

newLabels <- sapply(seqNames, function(x) taxonTable[which(taxonTable$sample == x)[1], 'manualMatch'])
table(newLabels == '' | is.na(newLabels))

newLabels <- gsub('\\s+', '_', newLabels)
names(newLabels) <- NULL
seqNames <- newLabels
seqNames <- setdiff(seqNames, 'Rhineura_floridana')
tax <- seqNames


# are any taxa not in the concatenated alignment?
tax <- intersect(tax, alnTaxa)

# put together genomic constraint tree
con1 <- keep.tip(genomTree, intersect(genomTree$tip.label, tax))
	
outgroups1 <- setdiff(outgroupTaxa, tax)
	
# write files
#- separate supermatrix for family
#- constraint tree, if needed
familyDir <- '~/Dropbox/squamatePhylo/2019/topologyScripts/datingTree'
if (!dir.exists(familyDir)) {
	dir.create(familyDir)
}

# alignment
## use phyx pxrms program to remove all but the specified taxa, then read in and write out to adjust ordering

# write a file containing the taxa to keep, and then use that to subset the alignment
write(c(outgroups1, tax), file = '~/Dropbox/squamatePhylo/2019/topologyScripts/datingTree/keeptaxa.txt')
args <- c('-s', fullAln, '-f', '~/Dropbox/squamatePhylo/2019/topologyScripts/datingTree/keeptaxa.txt', '-c', '-o', paste0(familyDir, '/famConcat.phy'))
system2('pxrms', args)
qq <- readDNAStringSet(paste0(familyDir, '/famConcat.phy'))
qq <- qq[c(outgroups1, tax)]
for (j in 1:length(qq)) {
	write(paste0('>', names(qq)[j]), file = paste0(familyDir, '/famConcat.aln'), append = ifelse(j == 1, FALSE, TRUE))
	write(as.character(qq[[j]]), file = paste0(familyDir, '/famConcat.aln'), append = TRUE)
}
file.remove('~/Dropbox/squamatePhylo/2019/topologyScripts/datingTree/keeptaxa.txt')
file.remove(paste0(familyDir, '/famConcat.phy'))

# constraint tree if needed
write.tree(con1, paste0(familyDir, '/famConstraint.tre'))

args <- c('iqtree2', '-s', 'famConcat.aln', '-spp', 'concatenated.partitions.iqtree', '-m', 'GTR+G', '-nt', nCores)

# add constraint	
args <- c(args, '-g', 'famConstraint.tre')

# add nIter
args <- c(args, '--runs', nIter, '-pre', 'datingTree_default')
	
# runtime
runTime <- "5:00:00:00"

# modify PBS template
PBScommand <- args
	
PBScommand[1] <- '/N/u/ptitle/Karst/iqtree2'
tmpPBS <- pbsTemplate
tmpPBS[grep('PBS -m a', tmpPBS)] <- gsub('-m a', '-m abe', grep('PBS -m a', tmpPBS, value=TRUE))
tmpPBS[grep('jobName', tmpPBS)] <- gsub('jobName', families[i], grep('jobName', tmpPBS, value=TRUE))
tmpPBS[grep('totalTime', tmpPBS)] <- gsub('totalTime', runTime, grep('totalTime', tmpPBS, value=TRUE))
tmpPBS[grep('tier1Families', tmpPBS)] <- 'cd /N/dc2/scratch/ptitle/datingTree;'
tmpPBS[grep('IQTREEcall', tmpPBS)] <- paste0(paste(PBScommand, collapse = ' '), ';')
tmpPBS <- tmpPBS[- grep('familyPath', tmpPBS)]

write(tmpPBS, file = '~/Dropbox/squamatePhylo/2019/topologyScripts/datingTree/datingTree.pbs')


