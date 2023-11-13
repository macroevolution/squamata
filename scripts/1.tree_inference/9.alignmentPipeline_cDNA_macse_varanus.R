
# to run such that output and errors are saved to file:
# Rscript --no-save --no-restore --verbose ~/Documents/pyphlawd/alignmentPipeline_cDNA.R all > outputFile.Rout 2>&1

# Run SuperCRUNCH scripts, MACSE, MAFFT, TCOFFEE, QuickProbs2, FASTTREE, RAxML

# nthreads is only defined for RAXML, all other programs use available cores.

# -----------------------------------------------------
dir1 <- './fasta_June2020_trimmedFiltered'

setwd(dir1)

faFiles <- list.files(pattern='\\.fa(sta)?$')
faFiles <- faFiles[!grepl('rRNA', faFiles)]

wd <- getwd()

availableThreads <- parallel::detectCores()
availableThreads <- 20

# only process loci with fewer than this many samples
lengthThreshold <- 3000

message(availableThreads, ' threads used.')

message(length(faFiles), ' files found.')

# Steps to run
runMACSE <- F
runAlignment <- T
runTreeBuilding <- F

message('\n\n-----------------------------------------')
message('\nPreparing to run the following steps:')
if (runMACSE) message('\tMACSE initial alignment')
if (runAlignment) message('\tprotein alignment')
if (runTreeBuilding) message('\ttree inference')
message('\n-----------------------------------------\n\n')

# -----------------------------------------------------
library(Biostrings)
require(parallel)
require(pbapply)

# -----------------------------------------------------
# define operations

# alignment method: 
## options:
## 	* mafft FFT (mafft_fft)
##	* FAMSA (FAMSA)
## 	* MAFFT G-large-ins-1 (mafft_large)
## 	* TCOFFEE with MAFFT (tcoffee_mafft)
## 	* TCOFFEE with QuickProbs2 (tcoffee_quickprobs2)
## 	* QuickProbs2 (only)
alnMethod <- c('mafft_large', 'tcoffee_mafft', 'tcoffee_quickprobs2', 'FAMSA', 'QuickProbs2')
alnMethod <- 'mafft_large'

# tree inference method: FastTree or RAxML-NG or RAxML
#treeMethod <- 'RAxML-NG'
#treeMethod <- 'FastTree'
treeMethod <- 'IQTREE_fast'

# run step1
runStep1 <- FALSE

# remove intermediate files and directories?
removeIntermediateFiles <- FALSE

# trim alignments?
doTrimming <- FALSE

# -----------------------------------------------------
# set some custom paths

# system <- 'mac'
# system <- 'varanus'
# system <- 'carbonate'

if (system2('echo', args='$HOME', stdout=T) == '/Users/pascaltitle') {
	superCrunchDir <- '~/SuperCRUNCH/supercrunch-scripts/'
	pythonPath <- '/usr/local/bin/python'
	macsePath <- '~/Downloads/macse_v2.03.jar'
	tcoffeePath <- '/Users/pascaltitle/tcoffee/bin/t_coffee'
	FastTreePath <- 'FastTree'
	raxmlNG_Path <- 'raxml-ng'
	raxmlPath <- 'raxmlHPC-PTHREADS-AVX'
	pxs2phyPath <- 'pxs2phy'
	pxaa2cdnPath <- 'pxaa2cdn'
	pxtlatePath <- 'pxtlate'
	famsaPath <- '/Users/pascaltitle/famsa-1.2.5-osx'
	quickProbsPath <- '~/Downloads/quickprobs-2.06-osx'
	javaMem <- '-Xmx10g' # amount of memory to allocate to java macse
}

if (system2('echo', args='$HOME', stdout=T) == '/home/ptitle@ads.iu.edu') {
	superCrunchDir <- '/home/ptitle@ads.iu.edu/SuperCRUNCH/supercrunch-scripts/'
	pythonPath <- 'python'
	macsePath <- '/home/ptitle@ads.iu.edu/Downloads/macse_v2.03.jar'
	tcoffeePath <- '/home/ptitle@ads.iu.edu/tcoffee/bin/t_coffee'
	FastTreePath <- 'FastTree'
	raxmlNG_Path <- '/home/ptitle@ads.iu.edu/Downloads/raxml-ng_v0.9.0/raxml-ng'
	raxmlPath <- 'raxmlHPC-PTHREADS-AVX'
	pxs2phyPath <- ''
	pxaa2cdnPath <- 'pxaa2cdn'
	pxtlatePath <- 'pxtlate'
	famsaPath <- '/home/ptitle@ads.iu.edu/famsa-1.2.5-linux-static'
	quickProbsPath <- '/home/ptitle@ads.iu.edu/Downloads/quickprobs-2.06-linux'
	javaMem <- '-Xmx40g' # amount of memory to allocate to java macse
}

if (system2('echo', args='$HOME', stdout=T) == '/N/u/ptitle/Carbonate') {
	superCrunchDir <- '/N/u/ptitle/Carbonate/SuperCRUNCH/supercrunch-scripts/'
	pythonPath <- 'python'
	macsePath <- '/N/u/ptitle/Carbonate/bin/macse_v2.03.jar'
	tcoffeePath <- '/N/u/ptitle/Carbonate/tcoffee/bin/t_coffee'
	FastTreePath <- '/N/u/ptitle/Carbonate/bin/FastTree'
	raxmlNG_Path <- 'raxml-ng'
	raxmlPath <- 'raxmlHPC-PTHREADS'
	pxs2phyPath <- ''
	pxaa2cdnPath <- '/N/u/ptitle/Carbonate/bin/pxaa2cdn'
	pxtlatePath <- '/N/u/ptitle/Carbonate/bin/pxtlate'
	quickProbsPath <- '/N/u/ptitle/Carbonate/bin/quickprobs-2.06-linux'
	javaMem <- '-Xmx20g' # amount of memory to allocate to java macse
}

if (system2('echo', args='$HOME', stdout=T) == '/N/u/ptitle/Karst') {
	superCrunchDir <- '/N/u/ptitle/Carbonate/SuperCRUNCH/supercrunch-scripts/'
	pythonPath <- 'python'
	macsePath <- '/N/u/ptitle/Karst/bin/macse_v2.03.jar'
	tcoffeePath <- '/N/u/ptitle/Karst/tcoffee/bin/t_coffee'
	FastTreePath <- '/N/u/ptitle/Karst/bin/FastTree'
	raxmlNG_Path <- 'raxml-ng'
	raxmlPath <- 'raxmlHPC-PTHREADS'
	pxs2phyPath <- ''
	pxaa2cdnPath <- '/N/u/ptitle/Karst/bin/pxaa2cdn'
	pxtlatePath <- '/N/u/ptitle/Karst/bin/pxtlate'
	quickProbsPath <- '/N/u/ptitle/Karst/bin/quickprobs-2.06-linux'
	javaMem <- '-Xmx20g' # amount of memory to allocate to java macse
}


#  -----------------------------------------------------
# Check that programs are runnable

if (Sys.which(pythonPath) == '') {
	stop('python not found.')
}

if (!file.exists(macsePath)) {
	stop('MACSE not found.')
}

if (Sys.which(tcoffeePath) == '') {
	stop('T-Coffee not found.')
}

if (Sys.which(FastTreePath) == '' & treeMethod == 'FastTree') {
	stop('FastTree not found.')
}

if (Sys.which(raxmlNG_Path) == '' & treeMethod == 'RAxML-NG') {
	stop('RAxML-ng not found.')
}

if (Sys.which(pxaa2cdnPath) == '') {
	stop('pxaa2cdn not found.')
}

if (Sys.which(pxtlatePath) == '') {
	stop('pxtlate not found.')
}

# if (Sys.which(raxmlPath) == '' & treeMethod == 'RAxML') {
	# stop('RAxML not found.')
# }

# if (Sys.which(pxs2phyPath) == '' & treeMethod == 'RAxML') {
	# stop('phyx pxs2phy not found.')
# }


# SuperCRUNCH Coding Translation Test to determine reading frame, and to separate seqs into more/less reliable
runCodingTranslationTest <- function(fasta, dnaType) {
	
	currentDir <- getwd()
	dir.create('orthologyTemp')
	file.copy(fasta, paste0('orthologyTemp/', fasta))
	args <- c(paste0(superCrunchDir, 'Coding_Translation_Tests.py'), '-i', paste0(currentDir, '/orthologyTemp'), '-o', paste0(currentDir, '/orthologyTemp'), '--table', dnaType, '--quiet')
	system2(command = pythonPath, args = args, stdout = '', stderr = '')
	
	# this will create outputs in several folders:  Translation-All-Seqs, Translation-Failed-Seqs, Translation-Passed-Seqs
	# move those to the main directory and delete
	file.copy(list.files(paste0(currentDir, '/orthologyTemp/Translation-All-Seqs'), pattern='\\.fasta', full.names=TRUE), to = './')
	if (length(list.files(paste0(currentDir, '/orthologyTemp/Translation-Passed-Seqs'), pattern='\\.fasta')) > 0) {
		file.copy(list.files(paste0(currentDir, '/orthologyTemp/Translation-Passed-Seqs'), pattern='\\.fasta', full.names=TRUE), to = './')
	}
	if (length(list.files(paste0(currentDir, '/orthologyTemp/Translation-Failed-Seqs'), pattern='\\.fasta')) > 0) {
		file.copy(list.files(paste0(currentDir, '/orthologyTemp/Translation-Failed-Seqs'), pattern='\\.fasta', full.names=TRUE), to = './')
	}
	
	system(paste0('rm -rf ', currentDir, '/orthologyTemp'))
}



# This function will do an initial MACSE alignment
macseAlign <- function(faFile) {
	
	if (grepl(paste0(c('^COI', '^CYTB', '^ND1', '^ND2', '^ND4'), collapse='|'), faFile)) {
		dnaType <- 'vertmtdna'
	} else {
		dnaType <- 'standard'
	}
	if (grepl('12S|16S', faFile)) {
		dnaType <- 'noncoding'
	}
	
	tmpDir <- paste0(gsub('\\.fa(sta)?', '', faFile), '_temp')
	if (dir.exists(tmpDir)) {
		stop('dir ', tmpDir, ' already exists.')
	}
	
	dir.create(tmpDir)
	
	# copy fasta file into temp dir, rename to .fasta and remove any underscores in name
	workingFaFile <- gsub('_', '-', faFile)
	workingFaFile <- gsub('\\.fa(sta)?', '.fasta', workingFaFile)
	file.copy(faFile, to = paste0(tmpDir, '/', workingFaFile))

	workingDir <- paste0(wd, '/', tmpDir)
	
	gene <- gsub('\\.fa(sta)?', '', faFile)

	step2Dir <- paste0(workingDir, '/step2')
	dir.create(step2Dir)
	setwd(step2Dir)
	
	file.copy(paste0(workingDir, '/', workingFaFile), './seqs.fa')
	
	# run orthology test
	runCodingTranslationTest('seqs.fa', dnaType)

	# Some genes have been causing MACSE to crash, so we will avoid that by splitting those genes in two, and merging the results. 
	if (grepl('COI|CYTB|cmos|CAND1|^ND2|^ND4', gene, ignore.case = TRUE)) {
		
		# try splitting the gene in two
		if (file.exists('seqs_Passed.fasta')) {
			seq1 <- scan('seqs_Passed.fasta', what='character', sep = '\n', quiet=TRUE)
			pass <- as.list(seq1[seq(2, length(seq1), by = 2)])
			names(pass) <- seq1[seq(1, length(seq1), by = 2)]
			rm(seq1)
		}
		if (file.exists('seqs_Failed.fasta')) {
			seq1 <- scan('seqs_Failed.fasta', what='character', sep = '\n', quiet=TRUE)
			fail <- as.list(seq1[seq(2, length(seq1), by = 2)])
			names(fail) <- seq1[seq(1, length(seq1), by = 2)]
			rm(seq1)
		}
		seq1 <- scan('seqs_All.fasta', what='character', sep = '\n', quiet=TRUE)
		all <- as.list(seq1[seq(2, length(seq1), by = 2)])
		names(all) <- seq1[seq(1, length(seq1), by = 2)]
		rm(seq1)
	
		# write first half
		fn <- 'seqs_All_A.fasta'
		for (j in 1:ceiling(length(all)/2)) {
			write(names(all)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
			write(all[[j]], file = fn, append = TRUE)
		}
	
		if (file.exists('seqs_Passed.fasta')) {
			fn <- 'seqs_Passed_A.fasta'
			for (j in 1:ceiling(length(pass)/2)) {
				write(names(pass)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
				write(pass[[j]], file = fn, append = TRUE)
			}
		}
		
		if (file.exists('seqs_Failed.fasta')) {
			fn <- 'seqs_Failed_A.fasta'
			for (j in 1:ceiling(length(fail)/2)) {
				write(names(fail)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
				write(fail[[j]], file = fn, append = TRUE)
			}
		}
	
		# write second half
		fn <- 'seqs_All_B.fasta'
		for (j in (ceiling(length(all)/2) + 1): length(all)) {
			write(names(all)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
			write(all[[j]], file = fn, append = TRUE)
		}
	
		if (file.exists('seqs_Passed.fasta')) {
			fn <- 'seqs_Passed_B.fasta'
			for (j in (ceiling(length(pass)/2) + 1): length(pass)) {
				write(names(pass)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
				write(pass[[j]], file = fn, append = TRUE)
			}
		}
		
		if (file.exists('seqs_Failed.fasta')) {
			fn <- 'seqs_Failed_B.fasta'
			for (j in (ceiling(length(fail)/2) + 1): length(fail)) {
				write(names(fail)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
				write(fail[[j]], file = fn, append = TRUE)
			}
		}
	
		# faster implementation to just get frameshift identification
		args <- c('-jar', javaMem, macsePath, '-prog', 'alignSequences', '-gc_def', ifelse(dnaType == 'standard', 1, 2), '-seq', 'seqs_Passed_A.fasta', '-seq_lr', 'seqs_Failed_A.fasta', '-max_refine_iter', 0, '-out_AA', 'step2_AA_A.aln', '-out_NT', 'step2_NT_A.aln')
	
		if (!any(grepl('_Failed_A.fasta', list.files()))) {
			args <- args[-grep('-seq_lr|seqs_Failed_A.fasta', args)]
			args <- gsub('seqs_Passed_A.fasta', 'seqs_All_A.fasta', args)
		}
		
		if (!any(grepl('_Passed_A.fasta', list.files()))) {
			# all seqs failed orthology test. Provide params used for less reliable seqs: -fs_lr 15 -fs_lr_term 20 -stop_lr 40
			args <- args[-grep('-seq_lr|seqs_Failed_A.fasta', args)]
			args <- gsub('seqs_Passed_A.fasta', 'seqs_All_A.fasta', args)
			args <- c(args, '-fs', 15, '-fs_term', 20, '-stop', 40)
		}
	
		system2(command = 'java', args = args)
		
		args <- c('-jar', '-Xmx20g', macsePath, '-prog', 'alignSequences', '-gc_def', ifelse(dnaType == 'standard', 1, 2), '-seq', 'seqs_Passed_B.fasta', '-seq_lr', 'seqs_Failed_B.fasta', '-max_refine_iter', 0, '-out_AA', 'step2_AA_B.aln', '-out_NT', 'step2_NT_B.aln')
	
		if (!any(grepl('_Failed_B.fasta', list.files()))) {
			args <- args[-grep('-seq_lr|seqs_Failed_B.fasta', args)]
			args <- gsub('seqs_Passed_B.fasta', 'seqs_All_B.fasta', args)
		}
		
		if (!any(grepl('_Passed_B.fasta', list.files()))) {
			# all seqs failed orthology test. Provide params used for less reliable seqs: -fs_lr 15 -fs_lr_term 20 -stop_lr 40
			args <- args[-grep('-seq_lr|seqs_Failed_B.fasta', args)]
			args <- gsub('seqs_Passed_B.fasta', 'seqs_All_B.fasta', args)
			args <- c(args, '-fs', 15, '-fs_term', 20, '-stop', 40)
		}
		
		system2(command = 'java', args = args)
		
		# combine results
		seq1 <- scan('step2_AA_A.aln', what='character', sep = '\n', quiet=TRUE)
		resA <- as.list(seq1[seq(2, length(seq1), by = 2)])
		names(resA) <- seq1[seq(1, length(seq1), by = 2)]
		rm(seq1)
		seq1 <- scan('step2_AA_B.aln', what='character', sep = '\n', quiet=TRUE)
		resB <- as.list(seq1[seq(2, length(seq1), by = 2)])
		names(resB) <- seq1[seq(1, length(seq1), by = 2)]
		rm(seq1)
		
		combined <- c(resA, resB)
		fn <- 'step2_AA.aln'
		for (j in 1:length(combined)) {
			write(names(combined)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
			write(combined[[j]], file = fn, append = TRUE)
		}
	
		seq1 <- scan('step2_NT_A.aln', what='character', sep = '\n', quiet=TRUE)
		resA <- as.list(seq1[seq(2, length(seq1), by = 2)])
		names(resA) <- seq1[seq(1, length(seq1), by = 2)]
		rm(seq1)
		seq1 <- scan('step2_NT_B.aln', what='character', sep = '\n', quiet=TRUE)
		resB <- as.list(seq1[seq(2, length(seq1), by = 2)])
		names(resB) <- seq1[seq(1, length(seq1), by = 2)]
		rm(seq1)
		
		combined <- c(resA, resB)
		fn <- 'step2_NT.aln'
		for (j in 1:length(combined)) {
			write(names(combined)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
			write(combined[[j]], file = fn, append = TRUE)
		}
	
	
	} else {
		
		# full alignment implementation -- unnecessary given that we are aiming to align amino acids later
		args <- c('-jar', javaMem, macsePath, '-prog', 'alignSequences', '-gc_def', ifelse(dnaType == 'standard', 1, 2), '-seq', 'seqs_Passed.fasta', '-seq_lr', 'seqs_Failed.fasta', '-out_AA', 'step2_AA.aln', '-out_NT', 'step2_NT.aln')
		
		# faster implementation to just get frameshift identification
		args <- c('-jar', javaMem, macsePath, '-prog', 'alignSequences', '-gc_def', ifelse(dnaType == 'standard', 1, 2), '-seq', 'seqs_Passed.fasta', '-seq_lr', 'seqs_Failed.fasta', '-max_refine_iter', 0, '-out_AA', 'step2_AA.aln', '-out_NT', 'step2_NT.aln')
		
		
		if (!any(grepl('_Failed.fasta', list.files()))) {
			args <- args[-grep('-seq_lr|seqs_Failed.fasta', args)]
			args <- gsub('seqs_Passed.fasta', 'seqs_All.fasta', args)
		}
		
		if (!any(grepl('_Passed.fasta', list.files()))) {
			# all seqs failed orthology test. Provide params used for less reliable seqs: -fs_lr 15 -fs_lr_term 20 -stop_lr 40
			args <- args[-grep('-seq_lr|seqs_Failed.fasta', args)]
			args <- gsub('seqs_Passed.fasta', 'seqs_All.fasta', args)
			args <- c(args, '-fs', 15, '-fs_term', 20, '-stop', 40)
		}
		
		system2(command = 'java', args = args, stdout = '', stderr = '')
	}

	#message('\n\n-----------------------------------------------')
	#message('STEP 3: Unalign MACSE aln and convert to AA')
	#message('-----------------------------------------------\n')
	
	# Replace ! with N
	# remove -
	# Replace * (stop codon) with X
	
	# Generate 2 files
	## -MACSE-generated AA alignment: replace ! with N and ungap.
	## -MACSE-generated NT alignment: replace ! with N and ungap.
	
	# Read in NT alignment from MACSE, replace ! with N, (replace * with NNN and delete gaps (not needed in NT alignment)).
	step2aln1 <- scan('step2_NT.aln', what='character', sep = '\n', quiet=TRUE)
	step2aln <- as.list(step2aln1[seq(2, length(step2aln1), by = 2)])
	names(step2aln) <- step2aln1[seq(1, length(step2aln1), by = 2)]
	rm(step2aln1)
	# ungap
	step2aln <- lapply(step2aln, function(x) gsub('\\-', '', x))
	step2aln <- lapply(step2aln, function(x) gsub('\\!', 'N', x))
	newfile <- 'step2_NT_ungapped.aln'
	for (i in 1:length(step2aln)) {
		write(names(step2aln)[i], file = newfile, append = ifelse(i == 1, FALSE, TRUE))
		write(step2aln[[i]], file = newfile, append = TRUE)
	}
	
	# translate ungapped NT alignment to AA
	# phyx NT -> AA tool
	args <- c('-s', 'step2_NT_ungapped.aln', '-t', ifelse(dnaType == 'standard', 'std', 'vmt'), '-o', 'translated_AA.fasta')
	system2('pxtlate', args = args)
	
	# Replace stop codon '*' with X in translated_AA.fasta
	step2aln1 <- scan('translated_AA.fasta', what='character', sep = '\n', quiet=TRUE)
	step2aln <- as.list(step2aln1[seq(2, length(step2aln1), by = 2)])
	names(step2aln) <- step2aln1[seq(1, length(step2aln1), by = 2)]
	rm(step2aln1)
	# define stop codons as something that won't trip up aligners
	step2aln <- lapply(step2aln, function(x) gsub('\\*', 'X', x))
	newfile <- 'translated_AA.fasta'
	for (i in 1:length(step2aln)) {
		write(names(step2aln)[i], file = newfile, append = ifelse(i == 1, FALSE, TRUE))
		write(step2aln[[i]], file = newfile, append = TRUE)
	}
	
	setwd(wd)
}


# sort fasta files by seq length
## helpful because we can then have the biggest ones go first, and diagnose issues that might come with fasta file size

seqLength <- numeric(length(faFiles))
for (i in 1:length(faFiles)) {
	seqLength[i] <- length(readDNAStringSet(faFiles[i]))
}

faFiles <- faFiles[order(seqLength, decreasing = TRUE)]
seqLength <- sort(seqLength, decreasing = TRUE)

# filter by n samples
faFiles <- faFiles[which(seqLength < lengthThreshold)]

if (length(faFiles) < length(seqLength)) {
	message(length(faFiles), ' being processed.')
}


###################
### RUN MACSE

if (runMACSE) {

	res <- pblapply(faFiles, function(x) macseAlign(x), cl = availableThreads)
}


####################
# RUN ALIGNMENTS

setwd(wd)

if (runAlignment) {
	
	message('\n\n------------------------------------')
	message('STEP 4: Align amino acid sequences')
	message('------------------------------------------\n')

	message('Running alignment with ', alnMethod, '.\n')

	for (i in 1:length(faFiles)) {
		
		message('\t', faFiles[i])
		
		gene <- gsub('\\.fa(sta)?', '', faFiles[i])
		
		tmpDir <- paste0(gsub('\\.fa(sta)?', '', faFiles[i]), '_temp')
		workingDir <- paste0(wd, '/', tmpDir)
	
		step4Dir <- paste0(workingDir, '/step4')
		if (!dir.exists(step4Dir)) dir.create(step4Dir)
		setwd(step4Dir)
		
		step2Dir <- paste0(workingDir, '/step2')
		file.copy(paste0(step2Dir, '/translated_AA.fasta'), './seqs.fa')	
	
		if (alnMethod == 'mafft_large') {
			args <- c('--large', '--localpair', '--thread', '-1', 'seqs.fa > seqs_LINS1.aln')
			system2(command = 'mafft', args = args, stdout = '', stderr = '')
			
			alignedAA <- 'seqs_LINS1.aln'
		
		} else if (alnMethod == 'FAMSA') {
			
			args <- c('seqs.fa', 'seqs_famsa.aln')
			system2(command = famsaPath, args = args)
			
			alignedAA <- 'seqs_famsa.aln'		
		
		} else if (alnMethod == 'QuickProbs2') {
	
			args <- c('seqs.fa', '-o', 'seqs_quickProbs2.aln')
			args <- c('seqs.fa', '-o', 'seqs_quickProbs2.aln --mem-limit 45000 -t 20')
			system2(command = quickProbsPath, args = args)
			
			alignedAA <- 'seqs_quickProbs2.aln'
	
		} else if (alnMethod == 'pasta') {
			args <- c('-i', 'seqs.fa', '-d', 'protein')
			system2(command ='run_pasta.py', args = args)
			
			alignedAA <- 'pastajob.marker001.seqs.fa.aln'
			file.rename(alignedAA, 'seqs_pasta.aln')
			alignedAA <- 'seqs_pasta.aln'
		}


		message('\n\n-----------------------------------------------')
		message('STEP 5: Derive NT from AA alignment')
		message('-----------------------------------------------\n')
	
		step5Dir <- paste0(workingDir, '/step5')
		if (!dir.exists(step5Dir)) dir.create(step5Dir)
		setwd(step5Dir)
		
		# The reportGapsAA2NT takes as input a FASTA file with unaligned nucleotide sequences and a FASTA file of aligned sequences that are the amino acid translations of the nucleotide ones. Each sequence should hence be present with the exact same name in both files and should be three time longer in the nucleotide file than in the amino acid file (ignoring gaps).
	
		# unaligned nucleotide sequence
		unalignedNT <- paste0(step2Dir, '/step2_NT_ungapped.aln')
		alignedAA <- paste0(step4Dir, '/', alignedAA)
		
		file.copy(unalignedNT, './')
		
		file.copy(alignedAA, './')
		
		#phyx program to convert AA alignment back to NT alignment, given AA alignment and unaligned NT fasta file.
		# pxaa2cdn -a AA_Alignment.fa -n Unaligned_Nucleotide.fa -o CDN_aln.fa
		args <- c('-a', basename(alignedAA), '-n', basename(unalignedNT), '-o', gsub('\\.aln$', '_NT.aln', basename(alignedAA)))
		system2(pxaa2cdnPath, args = args)

		alnFiles <- gsub('\\.aln$', '_NT.aln', basename(alignedAA))
	
		# Validation, quick check
		
		seq1 <- readDNAStringSet(unalignedNT)
		aln <- readDNAStringSet(alnFiles)
		
		length(seq1) == length(aln)
		message('Are there the same number of samples in orig fasta file and alignment?')
		message(length(seq1) == length(aln))
		
		seq1 <- lapply(seq1, function(x) as.character(x))
		aln <- lapply(aln, function(x) as.character(x))
		
		identical(sort(names(seq1)), sort(names(aln)))
		aln <- aln[names(seq1)]
		
		aln <- lapply(aln, function(x) gsub('\\-', '', x))
		
		check <- c()
		for (j in 1:length(seq1)) {
			check[j] <- identical(seq1[[j]], aln[[j]])
		}
		table(check) / length(check)
		message('How many samples have the same number of non-gap bases?')
		message((table(check) / length(check))*100, '%')


		# ------------------------------------
		## STEP 6
		# Move result files, and strip asterisks out of seq names (they were added by translation tests in step2)
		
		finalAln <- paste0(gene, '_', gsub('_NT\\.aln$', '.aln', gsub('seqs_', '', alnFiles)))
		
		aln <- readDNAStringSet(alnFiles)
		names(aln) <- gsub('\\s?\\*$', '', names(aln))
			
		for (j in 1:length(aln)) {
			write(paste0('>', names(aln)[j]), file = paste0(workingDir, '/', finalAln), append = ifelse(j == 1, FALSE, TRUE))
			write(as.character(aln[[j]]), file = paste0(workingDir, '/', finalAln), append = TRUE)
		}
		
		rm(aln)
		rm(finalAln)
		rm(unalignedNT)
		rm(alignedAA)
	}
}



	
	



