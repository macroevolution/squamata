# Take each locus alignment, and rename to taxon names

inputDir <- '~/Dropbox/Oz_Crown_Ages/finalAlignments_June2020'
outputDir <- '~/Dropbox/Oz_Crown_Ages/finalAlignments_June2020/renamed/'

require(Biostrings)

setwd(inputDir)

alnFiles <- list.files(pattern = '\\.aln$')

length(alnFiles)

tag <- '_onedirection.trimmedFiltered_LINS1_trimmed_outliersRemoved.aln|_onedirection.filtered_LINS1_trimmed_outliersRemoved.aln'

# read in table that allows us to interpret tip labels
taxonTableFile <- '~/Dropbox/Oz_Crown_Ages/fasta_March2020/metatable2.csv'
if (file.exists(taxonTableFile)) {
	taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
} else {
	taxonTableFile <- '~/squam2020/metatable2.csv'
	taxonTable <- read.csv(taxonTableFile, stringsAsFactors=FALSE)
}
head(taxonTable)

for (i in 1:length(alnFiles)) {
	
	gene <- gsub(tag, '', alnFiles[i])
	
	message('\t', gene)
	
	seq <- readDNAStringSet(alnFiles[i])
	
	subdat <- taxonTable[which(taxonTable$locus == gene), ]
	
	newLabels <- sapply(names(seq), function(x) subdat[which(subdat$accession == x), 'repdbTaxon'])
	newLabels <- gsub('\\s+', '_', newLabels)
	
	outfile <- paste0(outputDir, gsub('\\.aln$', '_renamed.aln', alnFiles[i]))
	
	for (j in 1:length(seq)) {
		write(paste0('>', newLabels[j]), file = outfile, append = ifelse(j == 1, FALSE, TRUE))
		write(as.character(seq[[j]]), file = outfile, append = TRUE)
	}
}

