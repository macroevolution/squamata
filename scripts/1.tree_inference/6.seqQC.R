#!/usr/bin/env Rscript

# Quality control checks on unaligned sequences:
# 	* count of ambiguous characters
# 	* count of stop codons and frameshifts as inferred with MACSE initial alignment

# -----------------------------------------------------
# deal with command line arguments
faFiles = commandArgs(trailingOnly=TRUE)

if (length(faFiles) == 0) {
  stop("usage: Rscript seqQC.R locus.fa(sta) OR all")
}

if (faFiles[1] == 'all') {
	faFiles <- list.files(pattern='\\.fa(sta)?$')
}

nthreads <- 0
# nthreads <- 18

javaMem <- '-Xmx20g'

message('Detected ', length(faFiles), ' seq files.\n')

deleteTempFiles <- TRUE

if (system2('uname', stdout=T) == 'Darwin') {
	macsePath <- '~/Downloads/macse_v2.03.jar'
	superCrunchDir <- '~/SuperCRUNCH/supercrunch-scripts/'
	pythonPath <- '/usr/local/bin/python'
} else if (system2('uname', stdout=T) == 'Linux') {
	macsePath <- '/home/ptitle@ads.iu.edu/Downloads/macse_v2.03.jar'
	superCrunchDir <- '/home/ptitle@ads.iu.edu/SuperCRUNCH/supercrunch-scripts/'
	pythonPath <- 'python'
}


require(parallel)
require(pbapply)


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


ambCodes <- c('R','Y','S','W','K','M','B','D','H','V')


# set up parallel processing
if (nthreads == 0) {
	cores <- detectCores()
	if (cores > length(faFiles)) {
		cores <- length(faFiles)
	}
} else {
	cores <- nthreads
}

message('\tUsing ', cores, ' cores.\n')


runQC <- function(faFile) {
	
	if (grepl(paste0(c('^COI', '^CYTB', '^ND1', '^ND2', '^ND4'), collapse='|'), faFile)) {
		dnaType <- 'vertmtdna'
	} else {
		dnaType <- 'standard'
	}
	if (grepl('12S|16S', faFile)) {
		dnaType <- 'noncoding'
	}
	
	tmpDir <- paste0(gsub('\\.fa(sta)?', '', faFile), '_QCtemp')
	if (dir.exists(tmpDir)) {
		stop('dir ', tmpDir, ' already exists.')
	}
	
	dir.create(tmpDir)
	
	# copy fasta file into temp dir, rename to .fasta and remove any underscores in name
	file.copy(faFile, to = paste0(tmpDir, '/', 'seqs.fa'))
	
	wd <- getwd()
	workingDir <- paste0(wd, '/', tmpDir)
	
	setwd(workingDir)
	
	# from unaligned seqs:
	# 	* number/percent of ambiguous characters

	seq1 <- scan('seqs.fa', what='character', sep = '\n', quiet=TRUE)
	seq <- as.list(seq1[seq(2, length(seq1), by = 2)])
	names(seq) <- seq1[seq(1, length(seq1), by = 2)]
	names(seq) <- gsub('^>', '', names(seq))
	names(seq) <- gsub('\\s?\\*$', '', names(seq))
	rm(seq1)

	seq <- lapply(seq, function(x) strsplit(x, '')[[1]])
	ambCount <- sapply(seq, function(x) sum(x %in% ambCodes))
	seqLength <- lengths(seq)
	ambPercent <- ambCount / seqLength

	
	if (dnaType != 'noncoding') {
		# run orthology test
		runCodingTranslationTest('seqs.fa', dnaType)


		if (grepl('CYTB|cmos|CAND1|^ND2|ND4', faFile)) {
			
			# try splitting the gene in two
			seq1 <- scan('seqs_Passed.fasta', what='character', sep = '\n', quiet=TRUE)
			pass <- as.list(seq1[seq(2, length(seq1), by = 2)])
			names(pass) <- seq1[seq(1, length(seq1), by = 2)]
			rm(seq1)
			seq1 <- scan('seqs_Failed.fasta', what='character', sep = '\n', quiet=TRUE)
			fail <- as.list(seq1[seq(2, length(seq1), by = 2)])
			names(fail) <- seq1[seq(1, length(seq1), by = 2)]
			rm(seq1)

			# write first half
			fn <- 'seqs_Passed_A.fasta'
			for (j in 1:ceiling(length(pass)/2)) {
				write(names(pass)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
				write(pass[[j]], file = fn, append = TRUE)
			}
			
			fn <- 'seqs_Failed_A.fasta'
			for (j in 1:ceiling(length(fail)/2)) {
				write(names(fail)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
				write(fail[[j]], file = fn, append = TRUE)
			}

			# write second half
			fn <- 'seqs_Passed_B.fasta'
			for (j in (ceiling(length(pass)/2) + 1): length(pass)) {
				write(names(pass)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
				write(pass[[j]], file = fn, append = TRUE)
			}
			
			fn <- 'seqs_Failed_B.fasta'
			for (j in (ceiling(length(fail)/2) + 1): length(fail)) {
				write(names(fail)[j], file = fn, append = ifelse(j == 1, FALSE, TRUE))
				write(fail[[j]], file = fn, append = TRUE)
			}

			# faster implementation to just get frameshift identification
			args <- c('-jar', javaMem, macsePath, '-prog', 'alignSequences', '-gc_def', ifelse(dnaType == 'standard', 1, 2), '-seq', 'seqs_Passed_A.fasta', '-seq_lr', 'seqs_Failed_A.fasta', '-max_refine_iter', 0, '-out_AA', 'step2_AA_A.aln', '-out_NT', 'step2_NT_A.aln')
			
			system2(command = 'java', args = args)
			
			args <- c('-jar', '-Xmx20g', macsePath, '-prog', 'alignSequences', '-gc_def', ifelse(dnaType == 'standard', 1, 2), '-seq', 'seqs_Passed_B.fasta', '-seq_lr', 'seqs_Failed_B.fasta', '-max_refine_iter', 0, '-out_AA', 'step2_AA_B.aln', '-out_NT', 'step2_NT_B.aln')
			
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
		
		} else {
		
			# faster implementation to just get frameshift identification
			args <- c('-jar', '-Xmx20g', macsePath, '-prog', 'alignSequences', '-gc_def', ifelse(dnaType == 'standard', 1, 2), '-seq', 'seqs_Passed.fasta', '-seq_lr', 'seqs_Failed.fasta', '-max_refine_iter', 0, '-out_AA', 'step2_AA.aln', '-out_NT', 'step2_NT.aln')
		
		
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
			
			system2(command = 'java', args = args)	
		}

		# from MACSE initial alignment
		# 	* number of contiguous frameshifts
		# 	* number of internal stop codons
		seq1 <- scan('step2_AA.aln', what='character', sep = '\n', quiet=TRUE)
		seq <- as.list(seq1[seq(2, length(seq1), by = 2)])
		names(seq) <- seq1[seq(1, length(seq1), by = 2)]
		names(seq) <- gsub('^>', '', names(seq))
		names(seq) <- gsub('\\s?\\*$', '', names(seq))
		rm(seq1)
		seq2 <- lapply(seq, function(x) strsplit(x, '')[[1]])
		nStopCodons <- sapply(seq2, function(x) sum(x == '*'))
		seq <- lapply(seq, function(x) gsub('\\!+', '!', x))
		seq2 <- lapply(seq, function(x) strsplit(x, '')[[1]])
		nFrameShifts <- sapply(seq2, function(x) sum(x == '!'))
		nStopCodons <- nStopCodons[names(seqLength)]
		nFrameShifts <- nFrameShifts[names(seqLength)]
	} else {
		nStopCodons <- NA
		nFrameShifts <- NA
	}
	
	
	tab <- matrix(nrow = length(seq), ncol = 6)
	colnames(tab) <- c('label', 'seqLength', 'nAmbig', 'percentAmbig', 'nStopCodons', 'nFrameShifts')
	tab <- as.data.frame(tab, stringsAsFactors=FALSE)
	tab[, 'label'] <- names(seqLength)
	tab[, 'seqLength'] <- seqLength
	tab[, 'nAmbig'] <- ambCount
	tab[, 'percentAmbig'] <- ambPercent
	tab[, 'nStopCodons'] <- nStopCodons
	tab[, 'nFrameShifts'] <- nFrameShifts
	
	fn <- paste0(wd, '/', gsub('\\.fa(sta)?', '_seqMetrics.csv', faFile))
	write.csv(tab, file = fn, row.names=FALSE)
	
	if (deleteTempFiles) {
		system(paste0('rm -rf ', workingDir))
	}
	
	setwd(wd)	
}

res <- pblapply(faFiles, function(x) runQC(x), cl = cores)



	
	

