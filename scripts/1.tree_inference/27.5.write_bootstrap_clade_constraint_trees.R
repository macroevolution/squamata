# Prep files for tier 2 inference
# create clade-specific constraint trees from IQTREE ultrafast bootstraps
# write unique topologies, but keep track of which topologies belong together

if (grepl('pascal', getwd())) {
  basedir <- '~/Dropbox/Oz_Crown_Ages/'
} else {
  basedir <- '~/Dropbox/squamate_tree/'
}

outdir <- paste0(basedir, "/", "tier2/")
wd = paste0(basedir, '/', 'data/all_bootstraps/')
setwd(wd)

require(ape)
require(Biostrings)

# collapse to unique topologies
# t = read.tree("all_IQtree.ufboot")
# t1 = lapply(t, unroot)
# t2 = unique.multiPhylo(t1)
# xx = attr(t2, "old.index")
# class(t2) <- "multiPhylo"
# write.tree(t2, "IQ.ufboot.unique")
# write.csv(xx, "IQ.ufboot.index")
#
# # which of the 1000 topologies does each unique topology belong to?
# bslink <- numeric(length(t1))
# for (i in 1:length(t1)) {
#  	topoCheck <- phangorn::RF.dist(t2, t1[[i]])
#  	topoMatch <- which(topoCheck == 0)
#  	bslink[i] <- topoMatch
#  	if (length(topoMatch) > 1) stop()
#  }
# sort(table(bslink), decreasing = TRUE) / length(t1)
# 20 of the 29 unique topologies appear more than once in the 1000 bootstraps
#
# # save small file with percent of posterior that each unique topology represents
# write.csv(as.matrix(table(bslink)/length(t1)), file = '~/Dropbox/Oz_Crown_Ages/tier2/percentTopologyInBootstrapPosterior.csv')
# 
# # read in unique topologies from genomic bootstrapping
# bs <- read.tree('IQ.ufboot.unique')
# 
# # rename the genomic trees
# # read in data table to replace names
# taxonTable <- read.csv('../../sqCL_taxa_repDBmatching.csv', stringsAsFactors=FALSE)
# # manual matching supercedes repDB auto matching, unless no manual match found
# taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'repDBtaxon']
# head(taxonTable)
# table(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch))
# # for the few individuals that don't have a good match
# # just match it to the existing species name
# taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'species']
# 
# # make sure all the genomic tips have a match!
# # just check one tree bc all the trees have same tips
# table(bs[[1]]$tip.label %in% taxonTable$sample)
# # the only missing tip is an outgroup sample, anyways
# bs[[1]]$tip.label[!bs[[1]]$tip.label %in% taxonTable$sample]
# 
# newbs = vector("list", length(bs))
# for (i in 1:length(bs)) {
#   genomTree = bs[[i]]
#   
#   # remove outgroups
#   genomOutgroups <- c('taeGut2', 'ISIS373002', 
#                       'UMFS-10956c', 'allMis1', 'chrPic1', 
#                       'galGal5', 'hg38', 'H20145a')
#   # genomOutgroups %in% genomTree$tip.label
#   genomOutgroups <- genomOutgroups[genomOutgroups %in% genomTree$tip.label]
#   genomTree <- root(genomTree, outgroup = genomOutgroups)
#   genomTree <- ladderize(genomTree)
#   genomTree <- drop.tip(genomTree, genomOutgroups)
#   
#   newLabels <- sapply(genomTree$tip.label, function(x) taxonTable[which(taxonTable$sample == x)[1], 'manualMatch'])
#   table(newLabels == '' | is.na(newLabels))
#   newLabels <- gsub('\\s+', '_', newLabels)
#   names(newLabels) <- NULL
#   genomTree$tip.label <- newLabels
#   
#   # remove a few bad individuals
#   genomTree <- drop.tip(genomTree, c('Bothrocophias_hyoprora',
#                                      'Corytophanes_hernandesii', 
#                                      'Lophognathus_gilberti',
#                                      'Rhineura_floridana'))

#		NEW SECTION corrects for duplicate tip names, and would lead to fewer unique topologies.
#			WOULD CONFLICT WITH CURRENT DATA - SO BE CAREFUL
		####////####
		# #   dups <- genomTree$tip.label[duplicated(genomTree$tip.label)]
		  # # for (j in 1:length(dups)) {
		  	# # n <- genomTree$edge[which(genomTree$edge[,2] == which(genomTree$tip.label == dups[j])[1]),1]
		  	# # message('\t', Ntip(extract.clade(genomTree, n)), ' tips')
		  # # }
		  # # They are all sister pairs. 
		  # genomTree <- drop.tip(genomTree, unlist(lapply(dups, function(x) which(genomTree$tip.label == x)[2:length(which(genomTree$tip.label == x))])))
		####////####
		
#   newbs[[i]] = genomTree
# }
# 
# class(newbs) <- "multiPhylo"
# write.tree(newbs, "IQ.ufboot.unique.renamed")

# Confirming whether or not these topologies are indeed unique
bsRF <- as.matrix(RF.dist(bs))
diag(bsRF) <- NA
newbsRF <- as.matrix(RF.dist(newbs))
diag(newbsRF) <- NA

which(bsRF == 0, arr.ind = T)
which(newbsRF == 0, arr.ind = T)


# read in best-estimate constraint trees 
# to get appropriate taxa per group
# but some families were run for both constraint & unconstraint
# whereas some were run for unconstraint only
bestTreeFiles <- list.files('../../tier1/fulltree/bestFamilyTrees/', 
                            pattern = '_con\\.tre$', full.names=TRUE)
bestTrees <- lapply(bestTreeFiles, function(x) read.tree(x))
names(bestTrees) <- gsub('_con\\.tre$', '', basename(bestTreeFiles))

# for each clade, subset the trees to those taxa, and store unique topologies
# conveniently, these best trees have the necessary outgroups
cladeUniqueBS <- vector('list', length(bestTrees))
cladeTaxaCt <- rep(NA, length(bestTrees))
names(cladeUniqueBS) <- names(bestTrees)
for (i in 1:length(bestTrees)) {
	message(i, ' - ', names(bestTrees)[i])
  cladeTaxaCt[i] = Ntip(bestTrees[[i]])
	cladeTaxa <- bestTrees[[i]]$tip.label
	
	all(cladeTaxa %in% newbs[[1]]$tip.label)
	setdiff(cladeTaxa, newbs[[1]]$tip.label)
	
	cladeTaxa <- intersect(cladeTaxa, newbs[[1]]$tip.label)
	
	cladeBS <- lapply(newbs, function(x) keep.tip(x, cladeTaxa))
	cladeBS <- lapply(cladeBS, unroot)
	
	class(cladeBS) <- 'multiPhylo'
	uniqueTopo <- unique.multiPhylo(cladeBS)
	cladeUniqueBS[[i]] <- uniqueTopo
	message('\t', length(uniqueTopo), ' topologies')
	
}

length(cladeUniqueBS)

# for each unique bootstrap tree, record which unique clade topology belongs with it
cladeMapping <- matrix(nrow = length(newbs), ncol = length(bestTrees))
colnames(cladeMapping) <- names(bestTrees)

for (i in 1:length(newbs)) {
	
	# for each clade j, determine which unique topology belongs in BS tree i
	for (j in 1:length(bestTrees)) {
		commonsp <- intersect(newbs[[i]]$tip.label, cladeUniqueBS[[j]][[1]]$tip.label)
		topoCheck <- sapply(cladeUniqueBS[[j]], function(x) all.equal.phylo(unroot(keep.tip(newbs[[i]], commonsp)), x, use.edge.length = FALSE))
		topoMatch <- which(topoCheck == TRUE)
		if (length(topoMatch) > 1) stop()
		cladeMapping[i, j] <- topoMatch
	}
}

# write unique clade topologies to file

# file structure
# create one folder for each bootsrap tree
# in each bootstrap folder, create a folder for each clade

setwd(outdir)
allCommands = "commands.txt"
for (i in 1:length(cladeUniqueBS)) {
	for (j in 1:length(cladeUniqueBS[[i]])) {
	  fam = names(cladeUniqueBS)[i]
		dirname <- paste0(outdir, '/families/', fam, '_uniqueTopo_', j)
		if (!dir.exists(dirname)) {
			dir.create(dirname)
		}
		fn <- paste0(dirname, '/constraint.tre')
		write.tree(cladeUniqueBS[[i]][[j]], fn)
		
		tree = cladeUniqueBS[[i]][[j]]
		tottax = cladeTaxaCt[i]
		tax = tree$tip.label
		nCores1 = 12
		nIter = 10
		
		# if constraint tree just has 6 tips
		# that means that it outgroups only
		if (length(tax) > 6) {
		  hasConstraint = TRUE
		} else {
		  hasConstraint = FALSE
		}
		
		# copy the alignment file over
		old.file = paste0('../tier1/families/', fam, "/famConcat.aln")
		new.file = paste0(dirname, "/famConcat.aln")
		file.copy(old.file, new.file)
		
		# print out iqtree commands
		args <- c('iqtree2', '-s', 'famConcat.aln', 
		          '-spp', '../concatenated.partitions.iqtree',
		          '-m', 'GTR+G', '-nt', nCores1)

		if (hasConstraint) {
		  args <- c(args, '-g', 'constraint.tre')
		}
		
		if (tottax < 500) {
		  args <- c(args, '--runs', nIter, '-pre', 
		            paste0(fam, '_default'))
		} else {
		  args <- c(args, '--seed')
		  args <- lapply(1:nIter, function(x) c(args, paste(rep(x,3), collapse = ''), '-pre', paste0(fam, '_default_', x)))
		}
		
		if (is.list(args)) {
		  write(paste("cd ", dirname), file = allCommands, append = TRUE)
		  for (xx in 1:length(args)) {
		    cmd = paste(args[[xx]], collapse = " ")
		    write(paste("cd ", dirname), file = allCommands, append = TRUE)
		    write(cmd, file = allCommands, append = TRUE)
		  }
		  } else {
		  write(paste("cd ", dirname), file = allCommands, append = TRUE)
		  cmd = paste(args, collapse = " ")
		  write(cmd, file = allCommands, append = TRUE)
		}
	}		
}


# need to create a phylip file -- should just copy over the file
# need to write output commands

# write clade mapping table to file
write.csv(cladeMapping, 'cladeMapping.csv', row.names = FALSE)


