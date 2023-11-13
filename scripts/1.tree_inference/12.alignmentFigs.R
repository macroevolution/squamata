# plot alignments next to tree

setwd('~/squam2020/bestGeneTrees')

clusters <- list.files(getwd(), pattern = '\\.tre$|\\.treefile$')
clusters <- setdiff(clusters, grep('ncbi\\.tre', clusters, value=TRUE))

clusters <- gsub('\\.tre$|\\.treefile$', '', clusters)

message('\t', length(clusters), ' trees detected.')

if (!dir.exists('../alignmentFigs')) {
	dir.create('../alignmentFigs')
}

outdir <- '../alignmentFigs'

library(ape)
library(Biostrings)

toniniTreeFile <- '~/Documents/pyphlawd/squam_shl_new_Consensus_9755.tre'
if (file.exists(toniniTreeFile)) {
	tonini <- read.tree(toniniTreeFile)
} else {
	toniniTreeFile <- '~/Downloads/Tonini_Squamata_Tree_Methods_Data/squam_shl_new_Consensus_9755.tre'
}
tonini <- read.tree(toniniTreeFile)
tonini$tip.label <- gsub('_', ' ', tonini$tip.label)

# read in table that allows us to interpret tip labels
metaTableFile <- '~/Dropbox/Oz_Crown_Ages/fasta_March2020/metatable2.csv'
if (file.exists(metaTableFile)) {
	taxonTable <- read.csv(metaTableFile, stringsAsFactors=FALSE)
} else {
	metaTableFile <- '~/squam2020/metatable2.csv'
	taxonTable <- read.csv(metaTableFile, stringsAsFactors=FALSE)
}
head(taxonTable)

for (i in 1:length(clusters)) {

	root <- root.phylo

	geneName <- clusters[i]

	message('\t', geneName)

	tr <- read.tree(paste0(geneName, c('.tre', '.treefile'))[paste0(geneName, c('.tre', '.treefile')) %in% list.files()])
	#al <- read_alignment(alignmentfile)
	al <- readDNAStringSet(paste0(geneName, '.aln'))

	if (grepl('\\s', names(al)[1])) {
		qq <- names(al)
		qq <- strsplit(qq, ' ')
		qq <- sapply(qq, function(x) x[1])
		names(al) <- qq
	}

	if (grepl('^>', names(al)[2]) & !grepl('^>', names(al)[1])) {
		names(al)[1] <- paste0('>', names(al)[1])
	}

	if (any(grepl('^>', names(al)))) {
		names(al) <- gsub('^>', '', names(al))
	}

	if (any(grepl('^>', tr$tip.label))) {
		tr$tip.label <- gsub('^>', '', tr$tip.label)
	}

	if (any(grepl('@', names(al)))) {
		ind <- which(grepl('@', names(al)) == TRUE)
		qq <- strsplit(names(al)[ind], '@')
		qq <- sapply(qq, function(x) x[2])
		names(al)[ind] <- qq
	}

	if (any(grepl('@', tr$tip.label))) {
		ind <- which(grepl('@', tr$tip.label))
		qq <- strsplit(tr$tip.label[ind], '@')
		qq <- sapply(qq, function(x) x[2])
		tr$tip.label[ind] <- qq
	}


	# test
	# testSet <- sample(tr$tip.label, 20)
	# tr <- keep.tip(tr, testSet)
	# al <- al[testSet]

	# replace accession numbers with taxon names
	spNames <- sapply(tr$tip.label, function(x) taxonTable[which(taxonTable$accession == x)[1], 'repdbTaxon'])
	missingTreeTips <- which(is.na(spNames))
	tr$tip.label <- as.character(spNames)	

	spNames <- sapply(names(al), function(x) taxonTable[which(taxonTable$accession == x)[1], 'repdbTaxon'])
	missingAlnNames <- which(is.na(spNames))
	names(al) <- spNames

	# names(al) <- gsub('^_R_', '', names(al))
	# spNames <- sapply(names(al), function(x) taxonTable[which(taxonTable$accession == x)[1], 'taxon'])
	# missingAlnNames <- which(is.na(spNames))
	# names(al) <- spNames

	if (length(missingTreeTips) > 0 & length(missingAlnNames) > 0) {
		if (identical(sort(names(missingTreeTips)), sort(names(missingAlnNames)))) {
			tr <- drop.tip(tr, missingTreeTips)
			al <- al[setdiff(1:length(al), missingAlnNames)]
		}
	}

	anyNA(tr$tip.label)
	anyNA(names(al))

	length(intersect(tr$tip.label, names(al))) == length(al)

	if (any(grepl('Sphenodon', tr$tip.label))) {
		tr <- ladderize(root(tr, grep('Sphenodon', tr$tip.label, value=TRUE)))
	} else {

		# subset tonini tree to as many of these taxa as possible
		toniniSub <- drop.tip(tonini, setdiff(tonini$tip.label, tr$tip.label))

		# what is the basal most taxon?
		xx <- which(toniniSub$edge[,1] == Ntip(toniniSub) + 1)
		xx <- toniniSub$edge[xx, 2]

		if (length(geiger::tips(toniniSub, xx[1])) < length(geiger::tips(toniniSub, xx[2]))) {
			outgroup <- geiger::tips(toniniSub, xx[1])
		} else {
			outgroup <- geiger::tips(toniniSub, xx[2])
		}
		
		# reroot
		if (length(outgroup) == 1) {
			tr <- ladderize(root(tr, outgroup))
		} else {
			tr <- ladderize(root(tr, node=getMRCA(tr, outgroup)))
		}
	}

	# tr <- ladderize(tr)

	al <- al[rev(na.omit(tr$tip.label[tr$edge[,2]]))]
	
	
	almat <- lapply(al, as.character)
	almat <- lapply(almat, strsplit, split='')
	almat <- lapply(almat, function(x) x[[1]])
	almat <- do.call('rbind', almat)
#	almat <- almat == '-'
	almat <- tolower(almat)
	
	# # colors	
	# replace bases with numbers so as to be able to use a color ramp
	col <- alignfigR::define_palette('dna')
	names(col) <- names(alignfigR::define_palette('dna'))
	names(col) <- tolower(names(col))
	for (j in 1:length(col)) {
		qq <- which(almat == names(col)[j])
		almat[qq] <- j
	}
	col[grep('grey', col)] <- 'gray95'
	mode(almat) <- 'numeric'

	

	outputfile <- paste0(outdir, '/', geneName, '.png')
	
	png(outputfile, width=4, height=10, units='in', res=300, type='cairo')
	
	layout(matrix(1:2, nrow=1, ncol=2), widths=c(1/3, 2/3))
	plot.phylo(tr, show.tip.label=F, no.margin=T, edge.width = 0.25)
	plot.phylo(tr, edge.color='white', show.tip.label=F, no.margin=T)
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	xvals <- seq(lastPP$x.lim[1], lastPP$x.lim[2], length.out=ncol(almat))
	yvals <- sort(lastPP$yy[1:Ntip(tr)])
	image(xvals, yvals, t(almat[nrow(almat):1,]), col = col, useRaster=T, add=TRUE)
	
	mtext(paste(Ntip(tr), 'tips x', length(al[[1]]), 'sites'), side = 3, adj = 0.8, line = -1, cex=0.25)
	

	dev.off()


}
