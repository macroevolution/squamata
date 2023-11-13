# Trim alignments at the extremities: Use MACSE trimAlignment to drop sites that are 90% empty or greater. Move from the flanks inwards until this condition is no longer true. 

# We will move through each gene directory, locate the alignment, trim it and then move it to a final alignments directory. 

newAlnDir <- '~/squam2020/fasta_June2020_alignments'
macsePath <- '/home/ptitle@ads.iu.edu/Downloads/macse_v2.03.jar'
# -----------------------------------------------------
dir1 <- '~/Dropbox/squamatePhylo/2019/topologyScripts/fasta_June2020_trimmedFiltered'
dir2 <- '~/squam2020/fasta_June2020_trimmedFiltered'
if (dir.exists(dir1)) {
	setwd(dir1)
} else if (dir.exists(dir2)) {
	setwd(dir2)
} else {
	stop('dirs not found.')
}

faFiles <- list.files(pattern='\\.fa(sta)?$')

faFiles <- faFiles[!grepl('rRNA', faFiles)]

alignmentTag <- 'LINS1\\.aln$'

wd <- getwd()

for (i in 1:length(faFiles)) {
	
	message('\t', faFiles[i])
	
	gene <- gsub('\\.fa(sta)?', '', faFiles[i])
	
	tmpDir <- paste0(gsub('\\.fa(sta)?', '', faFiles[i]), '_temp')
	workingDir <- paste0(wd, '/', tmpDir)

	setwd(workingDir)
	alnFile <- list.files(pattern = alignmentTag)
	
	outfile <- gsub('\\.aln', '_trimmed.aln', alnFile)
	outfile <- paste0(newAlnDir, '/', outfile)
	
	# java -jar macse.jar -prog trimAlignment -align AMBN_all_mafft_refined.fasta -min_percent_NT_at_ends XX
	args <- c('-jar', macsePath, '-prog', 'trimAlignment', '-align', alnFile, '-min_percent_NT_at_ends', 0.1, '-out_NT', outfile)
		
	system2(command ='java', args = args, stdout = TRUE, stderr = TRUE)

}

# rRNA 12S/16S
macsePath <- '~/Downloads/macse_v2.03.jar'
setwd('~/Dropbox/squamatePhylo/2019/topologyScripts/rRNA_split')

faFiles <- list.files(pattern = 'muscle')
faFiles <- faFiles[!grepl('trimmed\\.aln$', faFiles)]

for (i in 1:length(faFiles)) {
	
	message('\t', faFiles[i])
	
	alnFile <- faFiles[i]
	
	outfile <- gsub('\\.aln', '_trimmed.aln', alnFile)
	
	# java -jar macse.jar -prog trimAlignment -align AMBN_all_mafft_refined.fasta -min_percent_NT_at_ends XX
	args <- c('-jar', macsePath, '-prog', 'trimAlignment', '-align', alnFile, '-min_percent_NT_at_ends', 0.1, '-out_NT', outfile)
		
	system2(command ='java', args = args, stdout = TRUE, stderr = TRUE)	
}

