# Post-processing of the fasta files sometimes leads to a sample being split into 2 genes (generally, neighboring mt loci). In those cases, the sample now appears in the new gene fasta, with the gene name appended to the accession. Here, we will seek those out, and update the metadata table accordingly. 

# We should be able to find these by identifying accessions that are not currently listed in the metadata table.

fastaDir <- '~/Dropbox/Oz_Crown_Ages/fasta_March2020'

require(Biostrings)

setwd(fastaDir)
dat <- read.csv('metatable.csv', stringsAsFactors=FALSE)
dat$postProc <- NA

newMetadataFile <- 'metatable2.csv'

# table for records that are being added (because gene was split)
newdat <- matrix(nrow = 0, ncol = ncol(dat))
colnames(newdat) <- colnames(dat)

faFiles <- c(list.files(pattern = 'onedirection.trimmed\\.fa$'), list.files(pattern = 'rRNA_\\d\\dS_onedirection.fa$'))
faFiles <- faFiles[- grep('RAG1_onedirection', faFiles)]
faFiles <- sort(faFiles)

length(faFiles)

for (i in 1:length(faFiles)) {
	
	gene <- gsub('_onedirection.trimmed.fa|_onedirection.fa', '', faFiles[i])
	message('\t', i, '\t', gene)
	seq <- readDNAStringSet(faFiles[i])
	
	# which seq labels no longer appear in the metadata table?
	subdat <- dat[which(dat$locus == gene), ]
	diff <- setdiff(names(seq), subdat$accession)
	diffacc <- gsub(paste0('_', gene, '$'), '', diff)
	# if (length(setdiff(diffacc, dat$accession)) > 0) stop()
	# diffacc[which(lengths(sapply(diffacc, function(x) grep(x, dat$accession))) == 0)]
	
	if (length(diff) > 0) {
		# if accession is already associated with focal locus, then the sample was trimmed, so update the accession label.
		# if accession is not associated with focal locus, then it must have been split off of another gene. Add record to this locus.
		 
		for (j in 1:length(diff)) {
			
			if (diffacc[j] %in% subdat$accession) {
				# accession was already part of gene, but now has an added gene to the label
				# presumably, this happened because some of the seq belongs to another gene
				ind <- which(dat$locus == gene & dat$accession == diffacc[j])
				dat[ind, 'accession'] <- diff[j]
				dat[ind, 'postProc'] <- 'name updated'
	
			} else {
				# accession is not in this locus currently. It must have been split off from another locus, and is being added here. 
				if (length(grep(diffacc[j], dat$accession)) > 0) {
					ind <- grep(diffacc[j], dat$accession)
					# accession might show up for multiple loci if this accession was already split up, and that's ok.
					newRec <- dat[ind[1], ]
					newRec[, c('locus','accession','seqLength')] <- NA
					newRec$locus <- gene
					newRec$accession <- diff[j]
					newRec$postProc <- 'added to this locus due to gene splitting'
					newdat <- rbind(newdat, newRec)
				} else {
					# accession is not showing up anywhere -- doesn't seem to happen. good.
					stop()		
				}	
			}
		}
	}
	
	# are there any accessions that were in the fasta file, but that have been removed?
	# This could be due to gene being a poor match, or to only a small seq length remaining post trim
	missing <- setdiff(subdat$accession, union(names(seq), gsub(paste0('_', gene, '$'), '', names(seq))))
	
	# mark these as to be ignored. they were either poor quality, or they now belong to other loci.
	if (length(missing) > 0) {
		for (j in 1:length(missing)) {
			ind <- which(dat$locus == gene & dat$accession == missing[j])
			if (length(ind) != 1) stop()
			dat[ind, 'postProc'] <- 'removed'
		}
	}
}

# All RAG1 seqs should be removed too (RAG1 was not touched as it is not a separate fasta file)
dat[which(dat$locus == 'RAG1'), 'postProc'] <- 'removed'
table(dat[which(dat$locus == 'RAG1'), 'postProc'], useNA = 'always')

# add new records in
dat2 <- rbind.data.frame(dat, newdat)


# now check that all seq labels are accounted for
allAcc <- vector('list', length(faFiles))
for (i in 1:length(faFiles)) {
	
	gene <- gsub('_onedirection.trimmed.fa|_onedirection.fa', '', faFiles[i])
	seq <- readDNAStringSet(faFiles[i])
	
	subdat <- dat2[which(dat2$locus == gene), ]
	if (!all(names(seq) %in% subdat$accession)) stop()

	allAcc[[i]] <- names(seq)
	names(allAcc)[i] <- gene
}

sum(lengths(allAcc))
table(dat2$postProc != 'removed' | is.na(dat2$postProc)) # checks out!

table(dat2$postProc, useNA = 'always')
dat2 <- dat2[ - which(dat2$postProc == 'removed'),]
nrow(dat2)

for (i in 1:length(faFiles)) {
	
	gene <- gsub('_onedirection.trimmed.fa|_onedirection.fa', '', faFiles[i])
	seq <- readDNAStringSet(faFiles[i])
	
	subdat <- dat2[which(dat2$locus == gene), ]
	if (!all(names(seq) %in% subdat$accession)) stop()
}

# update seqLength field
for (i in 1:length(faFiles)) {
	
	gene <- gsub('_onedirection.trimmed.fa|_onedirection.fa', '', faFiles[i])
	message('\t', i, '\t', gene)
	seq <- readDNAStringSet(faFiles[i])
	
	for (j in 1:length(seq)) {
		ind <- which(dat2$locus == gene & dat2$accession == names(seq)[j])
		if (length(ind) != 1) stop()
		dat2[ind, 'seqLength'] <- nchar(seq[[j]])
	}
}


dat2 <- dat2[with(dat2, order(family, genus, repdbTaxon, locus)), ]
rownames(dat2) <- NULL
head(dat2)

write.csv(dat2, newMetadataFile, row.names = FALSE)