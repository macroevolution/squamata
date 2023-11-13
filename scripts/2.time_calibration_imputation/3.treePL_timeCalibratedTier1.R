# # # Extract divergence times from MCMCtree analysis and format for treePL.
## Time calibrate full tree with treePL, by breaking down tree into a couple of subclades.

library(ape)
library(castor)
library(MCMCtreeR)

treeFile <- '~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/tier1/fulltree/bestTrees/fulltree_default_con_1_raxmlOpt.raxml.final.tre'

maindir <- '~/Dropbox/squamatePhylo/2019/treePL_fulltree/resubmission'
outdir <- paste0(maindir, '/crossValidation')

# iteration A has the best mean likelihood and best overall convergence properties of lnL and divergence times
backboneFile <- '~/Dropbox/squamatePhylo/2019/part10/withOutgroups/ddbd_ARS3_June2023/mcmctree_AR/calibS3v2/mcmctree_useData2/iterA/FigTree.tre'
mcmcPostFile <- '~/Dropbox/squamatePhylo/2019/part10/withOutgroups/ddbd_ARS3_June2023/mcmctree_AR/calibS3v2/mcmctree_useData2/iterA/mcmc.txt'

fossilfile <- '~/Dropbox/squamatePhylo/2019/squamate_fossil_calibrations_v5_20230522_MEHJ.csv'

# getSpanningTips
#	returns a pair of tips that span a given node. 
#	if the node is terminal, includes "NA"
getSpanningTips <- function(phy, node){	
	if (node <= length(phy$tip.label)){
		return(c(phy$tip.label[node], 'NA'));
	}else{
		dnodes <- phy$edge[,2][phy$edge[,1] == node];

		while (dnodes[1] > length(phy$tip.label)){
			dnodes[1] <- phy$edge[,2][phy$edge[,1] == dnodes[1]][1];
		}
		while (dnodes[2] > length(phy$tip.label)){
			dnodes[2] <- phy$edge[,2][phy$edge[,1] == dnodes[2]][1];
		}		
		dset <- phy$tip.label[dnodes];
		return(dset);
	}	
}

tree <- read.tree(treeFile)
tree <- root(tree, 'Sphenodon_punctatus', resolve.root = TRUE)
tree <- read.tree(text = write.tree(ladderize(tree), file = ''))
tree <- drop.tip(tree, 'Sphenodon_punctatus')

backboneTree <- readMCMCtree(backboneFile)
mcmcPost <- data.table::fread(mcmcPostFile, data.table = FALSE)

# are there any tip name issues?
setdiff(backboneTree[[1]]$tip.label, tree$tip.label) # just outgroups, good!

# backbone topology should be identical in full tree
all.equal.phylo(keep.tip(backboneTree[[1]], intersect(tree$tip.label, backboneTree[[1]]$tip.label)), keep.tip(tree, intersect(tree$tip.label, backboneTree[[1]]$tip.label)), use.edge.length = FALSE)



fossildat <- read.csv(fossilfile, stringsAsFactors=FALSE)

fossildat <- fossildat[which(fossildat[,1] != ''),]
fossildat <- fossildat[which(fossildat$use != ''),]
fossildat <- fossildat[which(fossildat$use == 'yes' & fossildat$revision == ''),]

# subset table to only fossils for which spanning taxa are in tree
useFossils <- c()

for (i in 1:nrow(fossildat)) {
	spanningTaxa <- strsplit(fossildat[i, 'spanning.taxa'], split = '\\|')[[1]]
	if (all(spanningTaxa %in% tree$tip.label)) {
		useFossils <- c(useFossils, i)
	}
}

# not in tree, but otherwise flagged for use
fossilsNotUsed <- fossildat[setdiff(1:nrow(fossildat), useFossils), ]

fossildat <- fossildat[useFossils, ]

# match up nodes in mcmcPost to nodes in full tree
mcmcCols <- grep('t_n', colnames(mcmcPost), value = TRUE)
mcmcCorrespondingNodes <- c()
for (i in 1:length(mcmcCols)) {
	if (all(geiger::tips(backboneTree[[1]], as.numeric(gsub('t_n', '', mcmcCols[i]))) %in% tree$tip.label)) {
		mcmcCorrespondingNodes[i] <- getMRCA(tree, geiger::tips(backboneTree[[1]], as.numeric(gsub('t_n', '', mcmcCols[i]))))
	} else {
		mcmcCorrespondingNodes[i] <- NA
	}
}


# organize: get all node ages from backbone tree ready for full tree
## Here, we will use median node ages as fixed values for treePL

## also add additional fossils that could not be used with the genomic backbone
remainingFossils <- fossildat[which(fossildat$tree == 'fullTree'),]

# confirm
nrow(remainingFossils) == length(which(sapply(fossildat$spanning.taxa, function(x) all(strsplit(x, split = '\\|')[[1]] %in% backboneTree[[1]]$tip.label)) == FALSE))


calibdat <- as.data.frame(matrix(nrow = length(mcmcCols), ncol = 6))
colnames(calibdat) <- c('MinAge', 'MaxAge', 'fossil', 'fulltreeNode', 'backboneNode', 'spanningTaxa')

calibdat2 <- as.data.frame(matrix(nrow = nrow(remainingFossils), ncol = 6))
colnames(calibdat2) <- c('MinAge', 'MaxAge', 'fossil', 'fulltreeNode', 'backboneNode', 'spanningTaxa')


for (i in 1:length(mcmcCols)) {
	
	calibdat[i, 'backboneNode'] <- as.numeric(gsub('t_n', '', mcmcCols[i]))
	calibdat[i, 'fulltreeNode'] <- mcmcCorrespondingNodes[i]
	calibdat[i, 'fossil'] <- 'secondaryCalibration'
	if (is.na(mcmcCorrespondingNodes[i])) {
		calibdat[i, 'spanningTaxa'] <- NA
	} else {	
		calibdat[i, 'spanningTaxa'] <- paste0(getSpanningTips(tree, mcmcCorrespondingNodes[i]), collapse = ' ')
	}
	calibdat[i, c('MinAge', 'MaxAge')] <- median(mcmcPost[, mcmcCols[i]]) * 100
	
}

for (i in 1:nrow(remainingFossils)) {
	
	calibdat2[i, 'backboneNode'] <- NA
	calibdat2[i, 'fulltreeNode'] <- getMRCA(tree, strsplit(remainingFossils[i, 'spanning.taxa'], split = '\\|')[[1]])
	calibdat2[i, 'spanningTaxa'] <- gsub('\\|', ' ', remainingFossils[i, 'spanning.taxa'])
	calibdat2[i, 'fossil'] <- remainingFossils[i, 'name']
	calibdat2[i, 'MinAge'] <- remainingFossils[i, 'min.age']
	calibdat2[i, 'MaxAge'] <- remainingFossils[i, 'max.age']
	if (is.na(calibdat2[i, 'MaxAge'])) {
		calibdat2[i, 'MaxAge'] <- calibdat2[i, 'MinAge']
	}
}

calibdat$MinAge <- as.numeric(calibdat$MinAge)
calibdat$MaxAge <- as.numeric(calibdat$MaxAge)
calibdat2$MinAge <- as.numeric(calibdat2$MinAge)
calibdat2$MaxAge <- as.numeric(calibdat2$MaxAge)

# drop those calibrations that do not have representatives in the tree (here, due to outgroups)
calibdat <- calibdat[!is.na(calibdat$fulltreeNode), ]

anyDuplicated(c(calibdat$fulltreeNode, calibdat2$fulltreeNode))

# Create control file for TreePL

adOpts <- 1:5
optadOpts <- 1:5
optcvad <- 1:5

treePLopts <- expand.grid(list(adOpts, optadOpts, optcvad))


# Define subtrees that will be calibrated separately
## Pick 1 species of each family that is also in genomic dataset.
## these will be used to identify the appropriate nodes, and to create the parred down tree for assembly. 

subtreeTips <- list(

	gekko = c('Coleonyx_variegatus', 'Crenadactylus_naso', 'Nephrurus_levis', 'Lialis_burtonis', 'Hemidactylus_frenatus', 'Gonatodes_humeralis', 'Thecadactylus_solimoensis'),
	
	scinco = c('Plestiodon_skiltonianus', 'Zonosaurus_laticaudatus', 'Smaug_mossambicus', 'Xantusia_vigilis'),
	
	laterata = c('Bipes_biporus', 'Monopeltis_guentheri', 'Dryadosaura_nordestina', 'Ameiva_ameiva', 'Lacerta_viridis', 'Rhineura_floridana', 'Ptychoglossus_brevifrontalis', 'Blanus_strauchi', 'Cadea_blanoides', 'Trogonophis_wiegmanni'),
	
	iguania = c('Chamaeleo_calyptratus', 'Moloch_horridus', 'Uma_notata', 'Dipsosaurus_dorsalis', 'Anolis_punctatus', 'Diplolaemus_darwinii', 'Polychrus_liogaster', 'Hoplocercus_spinosus', 'Leiocephalus_macropus', 'Phymaturus_palluma', 'Oplurus_cyclurus', 'Basiliscus_basiliscus', 'Crotaphytus_collaris', 'Stenocercus_squarrosus'),
	
	anguids = c('Heloderma_suspectum', 'Ophiodes_striatus', 'Anniella_pulchra', 'Elgaria_kingii', 'Xenosaurus_platyceps', 'Shinisaurus_crocodilurus', 'Varanus_griseus', 'Lanthanotus_borneensis'),
	
	snakes = c('Acrochordus_granulatus', 'Anilius_scytale', 'Liotyphlops_ternetzii', 'Charina_bottae', 'Casarea_dussumieri', 'Contia_tenuis', 'Cylindrophis_ruffus', 'Micrurus_fulvius', 'Gerrhopilus_mirus', 'Fordonia_leucobalia', 'Atractaspis_fallax', 'Namibiana_labialis', 'Loxocemus_bicolor', 'Python_regius', 'Trachyboa_boulengeri', 'Anilios_waitii', 'Rhinophis_saffragamus', 'Agkistrodon_contortrix', 'Xenopeltis_unicolor', 'Anomochilus_leonardi', 'Pareas_nuchalis', 'Kladirostratus_acutus', 'Fimbrios_klossi', 'Xenophidion_schaeferi', 'Xenotyphlops_grandidieri')
	
)

all(unlist(subtreeTips) %in% tree$tip.label)

sapply(subtreeTips, function(x) Ntip(extract.clade(tree, getMRCA(tree, x))))


# make sure all subtrees have a calibrated crown node

for (i in 1:length(subtreeTips)) {
	
	cladeTree <- extract.clade(tree, getMRCA(tree, subtreeTips[[i]]))
		
	# is cladeTree root a calibrated node? It needs to be.
	message('\t', (Ntip(cladeTree) + 1) %in% sapply(calibdat$spanningTaxa, function(x) {
		if (all(strsplit(x, '\\s+')[[1]] %in% cladeTree$tip.label)) {
			getMRCA(cladeTree, strsplit(x, '\\s+')[[1]])
		} else {
			0
		}
	}))
}

# write these subtrees
setwd(outdir)
for (i in 1:length(subtreeTips)) {
    
    cladeTree <- extract.clade(tree, getMRCA(tree, subtreeTips[[i]]))
    cladeTreeName <- paste0(names(subtreeTips)[i], '.tre')
    write.tree(cladeTree, paste0('../', cladeTreeName))
}
    

## We will do a first pass where only smooth = 10 is attempted with the cross validation. This will allow us to identify which opt/optad/optcvad values lead to unreasonably large values, without running the full cross validation.

setwd(outdir)

for (i in 1:length(subtreeTips)) {
	
	cladeTree <- extract.clade(tree, getMRCA(tree, subtreeTips[[i]]))
	
	head(calibdat)
		
	inClade <- c()
	for (j in 1:nrow(calibdat)) {
		inClade[j] <- all(strsplit(calibdat[j, 'spanningTaxa'], ' ')[[1]] %in% cladeTree$tip.label)
	}
	
	cladeDir <- paste0(outdir, '/', names(subtreeTips)[i])		
	if (!dir.exists(cladeDir)) dir.create(cladeDir)
	
	setwd(cladeDir)
	
	cladeTreeName <- paste0(names(subtreeTips)[i], '.tre')
#	write.tree(cladeTree, cladeTreeName)
	
	for (j in 1:nrow(treePLopts)) {
		
		rundir <- paste0('treeplOpt_', paste0(treePLopts[j, ], collapse = ''))
		
		if (!dir.exists(rundir)) {
			dir.create(rundir)
		}
	
		setwd(rundir)
	
		treefile <- cladeTreeName
		nthreads <- 1
		numsites <- 72152
		
		treePLfile <- 'treePL_tier1_fulltree_initialPass.txt'
		outfileTag <- paste0(names(subtreeTips)[i], '_treePL')
		
		write('# treePL config file for time calibration of tier1 phylogeny.', file = treePLfile, append = FALSE)
		write(paste('# written', Sys.time()), file = treePLfile, append = TRUE)
		
		write(paste0('\n\ntreefile = ../../../', treefile), file = treePLfile, append = TRUE)
		write(paste0('nthreads = ', nthreads), file = treePLfile, append = TRUE)
		write(paste0('numsites = ', numsites), file = treePLfile, append = TRUE)
		
		write('\n###############################\n', file = treePLfile, append = TRUE)
		
		for (k in 1:nrow(calibdat)) {
			
			if (inClade[k]) {
			
				write(paste0('mrca = calibNode_', calibdat[k, 'fulltreeNode'], ' ', calibdat[k, 'spanningTaxa']), file = treePLfile, append = TRUE)
				write(paste0('min = calibNode_', calibdat[k, 'fulltreeNode'], ' ', round(calibdat[k, 'MinAge'], 4)), file = treePLfile, append = TRUE)
				write(paste0('max = calibNode_', calibdat[k, 'fulltreeNode'], ' ', round(calibdat[k, 'MaxAge'], 4)), file = treePLfile, append = TRUE)
				write('', file = treePLfile, append = TRUE)
			}
		}
		
		write('###############################', file = treePLfile, append = TRUE)
		# write('# Additional fossil calibrations\n', file = treePLfile, append = TRUE)
		
		# for (j in 1:nrow(calibdat2)) {
			
			# if (!calibdat2[j, 'fossil'] %in% c('Coluber cadurci', 'Naja romani', 'Chthonophis subterraneus')) {
				# # fossils that conflict will be left out.
				
				# write(paste0('# ', calibdat2[j, 'fossil']), file = treePLfile, append = TRUE)
				# write(paste0('mrca = calibNode_', calibdat2[j, 'fulltreeNode'], ' ', calibdat2[j, 'spanningTaxa']), file = treePLfile, append = TRUE)
				# write(paste0('min = calibNode_', calibdat2[j, 'fulltreeNode'], ' ', round(calibdat2[j, 'MinAge'], 4)), file = treePLfile, append = TRUE)
				# write(paste0('max = calibNode_', calibdat2[j, 'fulltreeNode'], ' ', round(calibdat2[j, 'MaxAge'], 4)), file = treePLfile, append = TRUE)
				# write('', file = treePLfile, append = TRUE)
			# }
		# }
		
		write('\n###############################\n', file = treePLfile, append = TRUE)
		
		# initial parameter block for verifying that different options for opt/optad/optcvad identify the same smooth value
		write('cv', file = treePLfile, append = TRUE)
		write('randomcv', file = treePLfile, append = TRUE)
#		write('verbose', file = treePLfile, append = TRUE)
		write('thorough', file = treePLfile, append = TRUE)
		write(paste0('opt = ', treePLopts[j, 1]), file = treePLfile, append = TRUE)
		write(paste0('optad = ', treePLopts[j, 2]), file = treePLfile, append = TRUE)
		write(paste0('optcvad = ', treePLopts[j, 3]), file = treePLfile, append = TRUE)
#		write('log_pen', file = treePLfile, append = TRUE)
		
		write('cvstart = 10', file = treePLfile, append = TRUE)
		write('cvstop = 10', file = treePLfile, append = TRUE)
			
		write('\n', file = treePLfile, append = TRUE)
		
		write(paste0('outfile = ', outfileTag, '.tre'), file = treePLfile, append = TRUE)
		write(paste0('cvoutfile = ', outfileTag, '.cvout'), file = treePLfile, append = TRUE)
	
		setwd(cladeDir)
	}
	
	setwd(outdir)
	
}

# create folder and files for tree that will act as scaffolding for subtrees. This tree will contain the tips used to define the subtrees, PLUS all tips found outside these subtrees. This way, each subtree can be grafted. 
setwd(outdir)

otherTaxa <- setdiff(tree$tip.label, unlist(lapply(subtreeTips, function(x) extract.clade(tree, getMRCA(tree, x))$tip.label)))
cladeTree <- keep.tip(tree, c(otherTaxa, unlist(subtreeTips)))

any(otherTaxa %in% backboneTree[[1]]$tip.label)

# We will not include secondary calibrations from within the subtrees, since those would be arbitrarily included based on taxon choice

# get taxa that are in scaffold tree and that define crown nodes
scaffoldNodeTaxa <- lapply(subtreeTips, function(x) {
	
	allpairs <- expand.grid(list(x, x))
	allpairs[,1] <- as.character(allpairs[,1])
	allpairs[,2] <- as.character(allpairs[,2])
#	allpairs <- lapply(1:nrow(allpairs), function(i) as.character(sort(allpairs[i, ])))
	allpairs <- lapply(1:nrow(allpairs), function(i) sort(as.character(allpairs[i, ])))
	allpairs <- allpairs[which(sapply(allpairs, function(x) x[1] == x[2]) == FALSE)]
	allpairs <- unique(allpairs)
	
	allpairs[[which(sapply(allpairs, function(x) getMRCA(tree, x)) == getMRCA(tree, x))[1]]]
})
	
scaffoldNodeTaxa <- unlist(c(scaffoldNodeTaxa, otherTaxa))

# determine what node numbers this tree would contain in the full tree
allpairs <- expand.grid(list(scaffoldNodeTaxa, scaffoldNodeTaxa))
allpairs[,1] <- as.character(allpairs[,1])
allpairs[,2] <- as.character(allpairs[,2])
#	allpairs <- lapply(1:nrow(allpairs), function(i) as.character(sort(allpairs[i, ])))
allpairs <- lapply(1:nrow(allpairs), function(i) sort(as.character(allpairs[i, ])))
allpairs <- allpairs[which(sapply(allpairs, function(x) x[1] == x[2]) == FALSE)]
allpairs <- unique(allpairs)
scaffoldNodes <- sapply(allpairs, function(x) getMRCA(tree, x))

inClade <- calibdat$fulltreeNode %in% scaffoldNodes
sum(inClade)	

plot(cladeTree, cex = 0.4)
calibratedScaffoldNodes <- unique(sapply(allpairs[which((scaffoldNodes %in% calibdat$fulltreeNode) == TRUE)], function(x) getMRCA(cladeTree, x)))
nodelabels(rep('', length(calibratedScaffoldNodes)), node = calibratedScaffoldNodes, frame = 'circle', cex = 0.3, bg = 'orange')
	
cladeDir <- paste0(outdir, '/', 'scaffoldTree')		
if (!dir.exists(cladeDir)) dir.create(cladeDir)
	
setwd(cladeDir)
	
cladeTreeName <- 'scaffold.tre'
# write.tree(cladeTree, cladeTreeName)

gsub('crossValidation', '', outdir)
write.tree(cladeTree, paste0(gsub('crossValidation', '', outdir), 'scaffold.tre'))
	
for (j in 1:nrow(treePLopts)) {
	
	rundir <- paste0('treeplOpt_', paste0(treePLopts[j, ], collapse = ''))
	
	if (!dir.exists(rundir)) {
		dir.create(rundir)
	}

	setwd(rundir)

	treefile <- cladeTreeName
	nthreads <- 1
	numsites <- 72152
	
	treePLfile <- 'treePL_tier1_fulltree_fullcv.txt'
	outfileTag <- paste0(cladeTreeName, '_treePL')
	
	write('# treePL config file for time calibration of tier1 phylogeny.', file = treePLfile, append = FALSE)
	write(paste('# written', Sys.time()), file = treePLfile, append = TRUE)
	
	write(paste0('\n\ntreefile = ../../../', treefile), file = treePLfile, append = TRUE)
	write(paste0('nthreads = ', nthreads), file = treePLfile, append = TRUE)
	write(paste0('numsites = ', numsites), file = treePLfile, append = TRUE)
	
	write('\n###############################\n', file = treePLfile, append = TRUE)
	
	for (k in 1:nrow(calibdat)) {
		
		if (inClade[k]) {
			
			spanningpair <- allpairs[[which(scaffoldNodes == calibdat[k, 'fulltreeNode'])[1]]]
			spanningpair <- paste(spanningpair, collapse = ' ')
		
			write(paste0('mrca = calibNode_', calibdat[k, 'fulltreeNode'], ' ', spanningpair), file = treePLfile, append = TRUE)
			write(paste0('min = calibNode_', calibdat[k, 'fulltreeNode'], ' ', round(calibdat[k, 'MinAge'], 4)), file = treePLfile, append = TRUE)
			write(paste0('max = calibNode_', calibdat[k, 'fulltreeNode'], ' ', round(calibdat[k, 'MaxAge'], 4)), file = treePLfile, append = TRUE)
			write('', file = treePLfile, append = TRUE)
		}
	}
	
	write('###############################', file = treePLfile, append = TRUE)	
	write('\n###############################\n', file = treePLfile, append = TRUE)
	
	# initial parameter block for verifying that different options for opt/optad/optcvad identify the same smooth value
	write('cv', file = treePLfile, append = TRUE)
	write('randomcv', file = treePLfile, append = TRUE)
	write('verbose', file = treePLfile, append = TRUE)
	write('thorough', file = treePLfile, append = TRUE)
	write(paste0('opt = ', treePLopts[j, 1]), file = treePLfile, append = TRUE)
	write(paste0('optad = ', treePLopts[j, 2]), file = treePLfile, append = TRUE)
	write(paste0('optcvad = ', treePLopts[j, 3]), file = treePLfile, append = TRUE)
	write('log_pen', file = treePLfile, append = TRUE)
	
	write('\n', file = treePLfile, append = TRUE)
	
	write(paste0('outfile = ', outfileTag, '.tre'), file = treePLfile, append = TRUE)
	write(paste0('cvoutfile = ', outfileTag, '.cvout'), file = treePLfile, append = TRUE)

	setwd(cladeDir)
}



# create master execute file for all folders
setwd(outdir)
write('', file = 'runAlltreePL_initial.txt', append = FALSE)
for (i in 1:length(subtreeTips)) {
    
    cladeDir <- paste0('./', names(subtreeTips)[i])

    # create master execute
    write(apply(treePLopts, 1, function(x) {
        paste0('cd ', cladeDir, '; ', 'cd treeplOpt_', paste0(x, collapse = ''), '; treePL treePL_tier1_fulltree_initialPass.txt > output_initialPass.txt;')
    }), file = 'runAlltreePL_initial.txt', append = TRUE)
}

# add scaffold tree
cladeDir <- paste0('./', 'scaffoldTree')	
write(apply(treePLopts, 1, function(x) {
    paste0('cd ', cladeDir, '; ', 'cd treeplOpt_', paste0(x, collapse = ''), '; treePL treePL_tier1_fulltree_fullcv.txt > output_fullcv.txt;')
}), file = 'runAlltreePL_initial.txt', append = TRUE)





## run with parallel < runAlltreePL_initial.txt --bar --joblog treepl_initial.log --timeout 300


# Now we will check on the initial runs and determine which ones returned reasonable chi squared values
## for those that are valid, write a new config file and include it in the execute file

# First write the results to a table (so that we don't lose information if overwriting directories)

# Now read all .cvout files and store information in table
initialCVmat <- as.data.frame(matrix(nrow = nrow(treePLopts), ncol = length(subtreeTips) + 1))
colnames(initialCVmat) <- c('params', names(subtreeTips))
initialCVmat[, 1] <- apply(treePLopts, 1, function(x) paste0(x, collapse = ''))

# create object for all smoothing/crossvalidation
crossValResults <- replicate(7, as.data.frame(matrix(nrow = nrow(treePLopts), ncol = 6)), simplify = FALSE)
names(crossValResults) <- c(names(subtreeTips), 'scaffoldTree')
for (i in 1:length(crossValResults)) {
	crossValResults[[i]][,1] <- apply(treePLopts, 1, function(x) paste0(x, collapse = ''))
}

for (i in 1:length(subtreeTips)) {
	
	cladeDir <- paste0(outdir, '/', names(subtreeTips)[i])		
	setwd(cladeDir)

	ll <- list()
	for (j in 1:nrow(treePLopts)) {
		
		rundir <- paste0('treeplOpt_', paste0(treePLopts[j, ], collapse = ''))
		setwd(rundir)
		ff <- list.files(pattern = '\\.cvout$')
		ff <- ff[!grepl('allcv', ff)]
		if (length(ff) > 0) {
			dat <- scan(ff, what = 'character', sep = '\n', quiet = TRUE)
			ll[[j]] <- setNames(as.numeric(gsub('(chisq: \\(.+\\)\\s)(.+)', '\\2', dat)), gsub('(chisq: \\(.+\\)\\s)(.+)', '\\1', dat))
		} else {
			ll[j] <- 0
		}	
		
		setwd('../')
	}
	
	names(ll) <- apply(treePLopts, 1, function(x) paste0(x, collapse = ''))
	ll <- ll[lengths(ll) > 0]
	ll <- ll[ll != 0]
	
	for (j in 1:length(ll)) {
		initialCVmat[which(initialCVmat[,1] == names(ll)[j]), names(subtreeTips)[i]] <- ll[j]
	}
}

# scaffold -- because the scaffold is small, all smoothing values were run

cladeDir <- paste0(outdir, '/scaffoldTree')		
setwd(cladeDir)

ll <- list()
for (j in 1:nrow(treePLopts)) {
	
	rundir <- paste0('treeplOpt_', paste0(treePLopts[j, ], collapse = ''))
	setwd(rundir)
	ff <- list.files(pattern = '\\.cvout$')
	ff <- ff[!grepl('allcv', ff)]
	if (length(ff) > 0) {
		dat <- scan(list.files(pattern = '\\.cvout$'), what = 'character', sep = '\n', quiet = TRUE)
		ll[[j]] <- setNames(as.numeric(gsub('(chisq: \\(.+\\)\\s)(.+)', '\\2', dat)), gsub('(chisq: \\(.+\\)\\s)(.+)', '\\1', dat))
	} else {
		ll[j] <- 0
	}	
	
	setwd('../')
}

names(ll) <- apply(treePLopts, 1, function(x) paste0(x, collapse = ''))

ll <- ll[lengths(ll) == 5]

for (j in 1:length(ll)) {
	crossValResults[['scaffoldTree']][which(crossValResults[['scaffoldTree']][,1] == names(ll)[j]), 2:ncol(crossValResults[['scaffoldTree']])] <- ll[[j]]
}
	
# write initial cross val to file
setwd(maindir)
write.csv(initialCVmat, 'initialCrossValResults.csv', row.names = FALSE)

	


# # hard code best smoother, based on manual examination
# topSmoothVal <- c()
# topSmoothVal['gekko'] <- "chisq: (0.1) "
# topSmoothVal['scinco'] <- "chisq: (0.1) "
# topSmoothVal['laterata'] <- "chisq: (1000) "
# topSmoothVal['iguania'] <- "chisq: (0.1) "
# topSmoothVal['anguids'] <- "chisq: (0.1) "
# topSmoothVal['snakes'] <- "chisq: (0.1) "
# topSmoothVal['scaffold'] <- "chisq: (1) "

# # hard coded based on previous analyses
# allbestParams <- c()
# allbestParams['gekko'] <- "555"
# allbestParams['scinco'] <- "553"
# allbestParams['laterata'] <- "351"
# allbestParams['iguania'] <- "551"
# allbestParams['anguids'] <- "335"
# allbestParams['snakes'] <- "335"
# allbestParams['scaffold'] <- "332"

useHardCoded <- FALSE


#####

cvReps <- 5


setwd(outdir)
write('', file = 'runAlltreePL_fullcv.txt', append = FALSE)

for (i in 1:length(subtreeTips)) {
	
	setwd(outdir)
	cladeDir <- paste0(outdir, '/', names(subtreeTips)[i])
	setwd(cladeDir)
	
	xx <- initialCVmat[, names(subtreeTips)[i]]
	names(xx) <- initialCVmat$params
	# We just ran cross validation for just smooth of 10
	## Which ones returned a reasonable chi squared value?
	valid <- which(xx < 1e+30)
		
	message('\t', names(subtreeTips)[i], '\t', 100 * length(valid) / length(xx), '% valid')
	
	cladeTreeName <- paste0(names(subtreeTips)[i], '.tre')
	
	# write new config file if value is reasonable
	for (j in 1:length(valid)) {
		
		rundir <- paste0('treeplOpt_', names(valid[j]))
		
		if (!dir.exists(rundir)) {
			dir.create(rundir)
		}
	
		setwd(rundir)
		
		# read in config file and modify
		configfile <- scan('treePL_tier1_fulltree_initialPass.txt', quiet = TRUE, sep = '\n', what = 'character')
	
		# remove constraints on cvstart/cvstop
		configfile <- configfile[ - grep('cvstart|cvstop', configfile)]
		
		configfile <- c(configfile, 'log_pen')
		
		# change tags
		configfile[grep('cvoutfile', configfile)] <- paste0('cvoutfile = ', names(subtreeTips)[i], '_treePL_allcv.cvout')
		
		for (k in 1:cvReps) {
	
			# write new version
			configfile1 <- configfile
	
			# change tags
			configfile1[grep('cvoutfile', configfile1)] <- gsub('allcv.cvout', paste0('allcv_rep', k, '.cvout'), configfile1[grep('cvoutfile', configfile1)])
			
			write(configfile1, file = paste0('treePL_tier1_fulltree_fullcv_rep', k, '.txt'))	

		}
				
		setwd(cladeDir)
	}
	
	# now write new execute file
	for (k in 1:cvReps) {
		write(apply(treePLopts[valid, ], 1, function(x) {
		    paste0('cd ./', basename(cladeDir), '; ', 'cd treeplOpt_', paste0(x, collapse = ''), '; treePL treePL_tier1_fulltree_fullcv_rep', k, '.txt > output_fullcv_rep', k, '.txt;')
		}), file = '../runAlltreePL_fullcv.txt', append = TRUE)
	}
	
	
	# write(apply(treePLopts[valid, ], 1, function(x) {
	# 	paste0('cd treeplOpt_', paste0(x, collapse = ''), '; treePL treePL_tier1_fulltree_fullcv.txt > output_fullcv.txt;')
	# }), file = 'runAlltreePL_fullcv.txt')		

}

# run an execute file for all scaffold cross validations
write(apply(treePLopts, 1, function(x) {
	    paste0('cd ./', 'scaffoldTree', '; ', 'cd treeplOpt_', paste0(x, collapse = ''), '; treePL treePL_tier1_fulltree_fullcv.txt > output_fullcv.txt;')
	}), file = 'runAlltreePL_scaffold_fullcv.txt', append = FALSE)



## Now run full cross validation only for those parameter combinations that returned reasonable values.
## --timeout 1200 = 20 minute limit
## run with parallel < runAlltreePL_fullcv.txt --bar --joblog treepl_fullcv.log --timeout 1200

## Read in reduced cross validation results

# Due to some weirdness with cross validation resulting in ties between very different values, when both high and low smoothing values have similar support, given the sizes of the trees, we will preferentially pick the smaller smoother.


##########################################################

finaldir <- gsub('crossValidation', 'finalruns', outdir)

bestSmoothTier1 <- numeric(length(subtreeTips)+1)
names(bestSmoothTier1) <- c(names(subtreeTips), 'scaffold')

    

writeFiles <- T
reps <- paste0('rep', 1:50)


for (i in 1:length(subtreeTips)) {
	
	cladeDir <- paste0(outdir, '/', names(subtreeTips)[i])		
	setwd(cladeDir)
	
	message('\t', names(subtreeTips)[i])

	if (!useHardCoded) {
		chisq <- vector('list',  nrow(treePLopts))
		for (j in 1:nrow(treePLopts)) {
			
			rundir <- paste0('treeplOpt_', paste0(treePLopts[j, ], collapse = ''))
			setwd(rundir)
			ff <- list.files(pattern = 'allcv_rep\\d+\\.cvout$')
			if (length(ff) > 0) {
				tmp <- list()
				for (k in 1:length(ff)) {
					dat <- scan(ff[k], what = 'character', sep = '\n', quiet = TRUE)
					tmp[[k]] <- setNames(as.numeric(gsub('(chisq: \\(.+\\)\\s)(.+)', '\\2', dat)), gsub('(chisq: \\(.+\\)\\s)(.+)', '\\1', dat))
				}
				chisq[[j]] <- tmp
			} else {
				chisq[j] <- NA
			}	
			
			setwd('../')
		}
		names(chisq) <- apply(treePLopts, 1, function(x) paste0(x, collapse = ''))
		
		
		# some may still have had unreasonable values. Leave those out.
		chisq <- chisq[!unlist(lapply(chisq, is.null), recursive = FALSE)]
		chisq <- chisq[which(sapply(chisq, function(x) !all(is.na(x))) == TRUE)]
		chisq <- chisq[which(lengths(chisq) == 5)]
		chisq <- chisq[sapply(lapply(chisq, function(x) lengths(x)), function(x) all(x == 5))]
				
		# read in output files to get likelihoods
		ll <-  vector('list', length(chisq))
		names(ll) <- names(chisq)
		for (j in 1:length(chisq)) {
			rundir <- paste0('treeplOpt_', names(chisq)[j])
			setwd(rundir)
			ff <- list.files(pattern = 'output_fullcv_rep\\d+.txt')
			ll[[j]] <- sapply(ff, function(x) {
				out <- scan(x, what = 'character', sep = '\n', quiet = TRUE)
				if (grepl('(after opt calc\\d?:\\s)(\\d+(\\.\\d+)?)', tail(out, 1))) {
					as.numeric(gsub('(after opt calc\\d?:\\s)(\\d+(\\.\\d+)?)', '\\2', tail(out, 1)))
				} else {
					NA
				}
			})
			
			setwd(cladeDir)
		}
		
		for (j in 1:length(chisq)) {
			names(chisq[[j]]) <- paste0(names(chisq)[j], '_rep', 1:length(chisq[[j]]))
		}

		for (j in 1:length(ll)) {
			names(ll[[j]]) <- paste0(names(ll)[j], '_rep', 1:length(ll[[j]]))
		}

		chisq <- unlist(chisq, recursive = FALSE)
		names(chisq) <- gsub('^\\d+\\.', '', names(chisq))

		# ll <- unlist(ll, recursive = FALSE)
		# names(ll) <- gsub('^\\d+\\.', '', names(ll))
		
		valid <- c()
		for (j in 1:length(chisq)) {
			if (all(unlist(chisq[[j]]) < 1e+30)) { # 1e+30
				valid <- c(valid, j)
			}
		}
		chisq <- chisq[valid]
		
		
		best <- c()
		for (j in 1:length(chisq)) {
			# best[[j]] <- sapply(chisq[[j]], function(x) names(x)[which.min(x)])
			best[j] <- names(chisq[[j]])[which.min(chisq[[j]])]
		}
		names(best) <- names(chisq)
		
		table(unlist(best))
		
		
		## Write new config file with smoothing value and best opt/optad/optcvad parameter values
		
		## which was best smoothing value? (value that was most frequently returned)
		bestSmooth <- names(which(table(unlist(best)) == max(table(unlist(best)))))

		if (length(bestSmooth) > 1) {
			# choose the more conservative, larger value
			stop()
			bestSmooth <- bestSmooth[which.max(as.numeric(gsub('(chisq:\\s\\((\\d+(\\.\\d+)?)\\))', '\\2', bestSmooth)))]
		}
		runsWithBestSmooth <- which(unlist(best) == bestSmooth)

		bestParams <- names(runsWithBestSmooth)
		bestParams <- gsub('_rep\\d+$', '', bestParams)
		bestParams <- unique(bestParams)

	} else {
		
		# use hard-coded value
		bestSmooth <- topSmoothVal[names(subtreeTips)[i]]
		bestParams <- allbestParams[names(subtreeTips)[i]]
	}	
	
	bestSmooth <- as.numeric(gsub('(chisq:\\s\\((\\d+(\\.\\d+)?)\\))', '\\2', bestSmooth))
	bestSmoothTier1[names(subtreeTips)[i]] <- bestSmooth
	
	if (writeFiles) {
	
		counter <- 1
		for (j in 1:length(bestParams)) {
			cladeDir <- paste0(outdir, '/', names(subtreeTips)[i])		
			setwd(cladeDir)
		    config <- scan(paste0('treeplOpt_', bestParams[j], '/treePL_tier1_fulltree_fullcv.txt'), what = 'character', sep = '\n', quiet = TRUE)
	    
		    # make modifications to config
		    ## add smoothing value, remove cross validation flags
		    config[grep('treefile', config)] <- gsub('(treefile\\s=\\s)(../../../)(\\w+\\.tre)', '\\1../../\\3', config[grep('treefile', config)])
		    config <- config[ - grep('^randomcv$|^cv$|^cvoutfile', config)]
		    config <- c(config, '', paste0('smooth = ', bestSmooth))
		    
	    		# write
	    		cladeDir <- paste0(finaldir, '/', names(subtreeTips)[i])		
	    		if (!dir.exists(cladeDir)) dir.create(cladeDir)
	    
			setwd(cladeDir)
	    
		    for (k in 1:length(reps)) {
		        config1 <- config
		        config1[grep('outfile', config1)] <- gsub('\\.tre', paste0('_rep', counter, '.tre'), grep('outfile', config1, value = TRUE))
		        write(config1, file = paste0(names(subtreeTips)[i], '_treePLfinal_rep', counter, '.txt'))
		        counter <- counter + 1
		    }
		}
	}
	    
	setwd(outdir)
		
}


### Process scaffold tree cross validation results and prep final run files

cladeDir <- paste0(outdir, '/', 'scaffoldTree')		
setwd(cladeDir)

if (!useHardCoded) {
	
	chisq <- vector('list',  nrow(treePLopts))
	for (j in 1:nrow(treePLopts)) {
		
		rundir <- paste0('treeplOpt_', paste0(treePLopts[j, ], collapse = ''))
		setwd(rundir)
		if (length(list.files(pattern = '\\.cvout$')) > 0) {
			dat <- scan(list.files(pattern = '\\.cvout$'), what = 'character', sep = '\n', quiet = TRUE)
			chisq[[j]] <- setNames(as.numeric(gsub('(chisq: \\(.+\\)\\s)(.+)', '\\2', dat)), gsub('(chisq: \\(.+\\)\\s)(.+)', '\\1', dat))
		} else {
			chisq[i] <- 0
		}	
		
		setwd('../')
	}
	names(chisq) <- apply(treePLopts, 1, function(x) paste0(x, collapse = ''))
	
	# some may still have had unreasonable values. Leave those out.
	chisq <- chisq[!unlist(lapply(chisq, is.null))]
	chisq <- chisq[which(lengths(chisq) == 5)]
	
	valid <- c()
	for (j in 1:length(chisq)) {
		if (all(chisq[[j]] < 1e+30)) {
			valid <- c(valid, j)
		}
	}
	chisq <- chisq[valid]
	
	# read in output files to get likelihoods - we will use the last printed likelihood, which is from a full optimization.
	ll <-  numeric(length(chisq))
	names(ll) <- names(chisq)
	for (j in 1:length(chisq)) {
		rundir <- paste0('treeplOpt_', names(chisq)[j])
		setwd(rundir)
		
		out <- scan('output_fullcv.txt', what = 'character', sep = '\n', quiet = TRUE)
		if (grepl('(after opt calc:\\s)(\\d+(\\.\\d+)?)', tail(out, 1))) {
			ll[j] <- as.numeric(gsub('(after opt calc:\\s)(\\d+(\\.\\d+)?)', '\\2', tail(out, 1)))
		} else {
			ll[j] <- NA
		}		
		setwd(cladeDir)
	}
	
	
	best <- c()
	for (j in 1:length(chisq)) {
		best[j] <- names(chisq[[j]])[which.min(chisq[[j]])]
	}
	
	table(best)
	
	## which was best smoothing value? (value that was most frequently returned)
	bestSmooth <- names(which(table(best) == max(table(best))))	
	
	if (length(bestSmooth) > 1) {
		# choose the more conservative, larger value
		bestSmooth <- bestSmooth[which.max(as.numeric(gsub('(chisq:\\s\\((\\d+(\\.\\d+)?)\\))', '\\2', bestSmooth)))]
	}
	runsWithBestSmooth <- which(best == bestSmooth)
	# llWithBestSmooth <- ll[runsWithBestSmooth]
	# sort(llWithBestSmooth, decreasing = TRUE)
	# bestParams <- names(which.min(llWithBestSmooth))

	bestParams <- names(chisq)[runsWithBestSmooth]
	bestParams <- unique(bestParams)
	
	bestSmooth <- as.numeric(gsub('(chisq:\\s\\((\\d+(\\.\\d+)?)\\))', '\\2', bestSmooth))	

} else {
	
	# use hard-coded value
	bestSmooth <- as.numeric(gsub('(chisq:\\s\\((\\d+(\\.\\d+)?)\\))', '\\2', topSmoothVal['scaffold']))
	bestParams <- allbestParams['scaffold']
}

## Write new config file with smoothing value and best opt/optad/optcvad parameter values

bestSmoothTier1['scaffold'] <- bestSmooth

counter <- 1
for (i in 1:length(bestParams)) {

	cladeDir <- paste0(outdir, '/', 'scaffoldTree')		
	setwd(cladeDir)

	config <- scan(paste0('treeplOpt_', bestParams[i], '/treePL_tier1_fulltree_fullcv.txt'), what = 'character', sep = '\n', quiet = TRUE)

	# make modifications to config
	## add smoothing value, remove cross validation flags
	config[grep('treefile', config)] <- gsub('(treefile\\s=\\s)(../../../)(\\w+\\.tre)', '\\1../../\\3', config[grep('treefile', config)])
	config <- config[ - grep('^randomcv$|^cv$|^cvoutfile', config)]
	config <- c(config, '', paste0('smooth = ', bestSmooth))
	
	# write
	cladeDir <- paste0(finaldir, '/scaffoldTree')		
	if (!dir.exists(cladeDir)) dir.create(cladeDir)
	
	setwd(cladeDir)
	
	for (j in 1:length(reps)) {
	    config1 <- config
	    config1[grep('outfile', config1)] <- paste0('outfile = scaffold_treePL_rep', counter, '.tre') 
	    write(config1, file = paste0('scaffoldTree_treePLfinal_rep', counter, '.txt'))
	    counter <- counter + 1
	}
}

# write execution file
setwd(finaldir)
f <- list.files(recursive = TRUE, pattern = '\\.txt$')
write(paste0('cd ', sapply(strsplit(f, '/'), function(x) x[1]), '; treePL ', sapply(strsplit(f, '/'), function(x) x[2]), ' > output_', sapply(strsplit(f, 'rep'), function(x) x[2]), ';'), file = 'runAlltreePL_finalruns.txt')

# run with parallel --bar --joblog treePL.log < runAlltreePL_finalruns.txt




##############################
## Read in final treePL runs and assemble.

setwd(finaldir)

par(mfrow = c(3,2))

finalsubtrees <- vector('list', length(subtreeTips))
names(finalsubtrees) <- names(subtreeTips)
corList <- meandiffList
for (i in 1:length(subtreeTips)) {
	setwd(paste0('./', names(subtreeTips)[i]))
	ff <- list.files(pattern = 'output_\\d+\\.txt$')
    ll <-  numeric(length(ff))
    names(ll) <- paste0('rep', gsub('(output_)(\\d+)\\.txt', '\\2', ff))
    for (j in 1:length(ff)) {
        out <- scan(ff[j], what = 'character', sep = '\n', quiet = TRUE)
        ll[j] <- as.numeric(gsub('(after opt calc:\\s)(\\d+(\\.\\d+)?)', '\\2', tail(out, 1)))
    }
    plot(sort(ll)); mtext(paste0(names(subtreeTips)[i], ' ', min(ll)))
    message('\t', names(subtreeTips)[i], '\t', min(ll))
    message('\t', names(subtreeTips)[i], '\t', length(which(ll == min(ll))), ' with best ll')
    treeFile <- grep('rep\\d+\\.tre$', list.files(full.names = TRUE), value = TRUE)
	treeFile <- grep(paste0(names(ll)[which.min(ll)], '\\.tre'), treeFile, value = TRUE)
    finalsubtrees[[i]] <- read.tree(treeFile)
    
    # how different are the different runs?
    treeFile <- grep('rep\\d+\\.tre$', list.files(full.names = TRUE), value = TRUE)
    tr <- lapply(names(sort(ll)[1:10]), function(x) read.tree(grep(paste0(x, '\\.tre'), treeFile, value = TRUE)))
    sapply(tr[2:length(tr)], function(x) all.equal.phylo(tr[[1]], x, use.edge.length = FALSE))
    cors <- c()
    for (j in 2:length(tr)) {
	   cors[j] <- cor(branching.times(tr[[1]]), branching.times(tr[[j]]))
	}
	corList[[i]] <- cors[!is.na(cors)]
	
    setwd(finaldir)    
}

sapply(corList, range, simplify = FALSE)

# scaffold
setwd(paste0('./scaffoldTree'))
ff <- list.files(pattern = 'output_\\d+\\.txt$')
ll <-  numeric(length(ff))
names(ll) <- paste0('rep', gsub('(output_)(\\d+)\\.txt', '\\2', ff))
for (j in 1:length(ff)) {
    out <- scan(ff[j], what = 'character', sep = '\n', quiet = TRUE)
    ll[j] <- as.numeric(gsub('(after opt calc:\\s)(\\d+(\\.\\d+)?)', '\\2', tail(out, 1)))
}
message('\tscaffold\t', length(which(ll == min(ll))), ' with best ll')
plot(sort(ll))
treeFile <- grep('rep\\d+\\.tre$', list.files(full.names = TRUE), value = TRUE)
treeFile <- grep(paste0(names(ll)[which.min(ll)], '\\.tre'), treeFile, value = TRUE)
finalScaffold <- read.tree(treeFile)

# add tag to more easily be able to remove later
finalScaffold$tip.label <- paste0(finalScaffold$tip.label, '_scaffold')

# assemble
finalFullTree <- finalScaffold

# subtreeTips lists the subset of species that represent each subtree. We will swap out that subset with the full, time-calibrated tree that was calibrated separately. Crown ages should match.
for (i in 1:length(subtreeTips)) {

    finalScaffoldCopy <- finalScaffold
    finalScaffoldCopy$tip.label <- gsub('_scaffold$', '', finalScaffoldCopy$tip.label)
    branching.times(finalScaffoldCopy)[as.character(getMRCA(finalScaffoldCopy, subtreeTips[[i]]))]
    max(branching.times(finalsubtrees[[i]]))
	# they match. good.
	
	
	# replace finalFullTree with one where the subtree has been added in
	finalFullTreeCopy <- finalFullTree
	finalFullTreeCopy$tip.label <- gsub('_scaffold$', '', finalFullTreeCopy$tip.label)
	bindNode <- getMRCA(finalFullTreeCopy, subtreeTips[[i]])
	finalFullTree <- bind.tree(finalFullTree, finalsubtrees[[i]], where = bindNode, position = 0)
	toDrop <- paste0(finalsubtrees[[i]]$tip.label, '_scaffold')[which(paste0(finalsubtrees[[i]]$tip.label, '_scaffold') %in% finalFullTree$tip.label)]
	finalFullTree <- drop.tip(finalFullTree, toDrop)
	
	# branching time should again match.
	branching.times(finalFullTree)[as.character(getMRCA(finalFullTree, subtreeTips[[i]]))]

	# should be ultrametric still!
	is.ultrametric(finalFullTree, option = 2)
	is.binary(finalFullTree)
}

length(unique(finalFullTree$tip.label))
Ntip(finalFullTree)

grep('_scaffold$', finalFullTree$tip.label, value = TRUE)
finalFullTree$tip.label <- gsub('_scaffold$', '', finalFullTree$tip.label)


# final tree is not ultrametric, but only because of numerical precision. We will fix this with phytools::force.ultrametric.
is.ultrametric(finalFullTree)
is.ultrametric(finalFullTree, option = 2)

finalFullTree <- ape::chronoMPL(finalFullTree)
finalFullTree2 <- phytools::force.ultrametric(finalFullTree, method = 'nnls') # for comparison

# final check
is.ultrametric(finalFullTree)
is.binary(finalFullTree)
is.rooted(finalFullTree)
length(unique(finalFullTree$tip.label)) == Ntip(finalFullTree)
all.equal.phylo(finalFullTree, tree, use.edge.length = FALSE)

# do dates still match MCMCtree time-calibrated tree?
for (i in 1:nrow(calibdat)) {
    tt <- strsplit(calibdat[i, 'spanningTaxa'], '\\s+')[[1]]
    fulltreeAge <- branching.times(finalFullTree)[as.character(getMRCA(finalFullTree, tt))]
    backboneAge <- calibdat[i, 'MinAge']
    if ((abs(fulltreeAge - backboneAge) < 1e-4) == FALSE) message('\t', i)
}

plot.phylo(finalFullTree, show.tip.label = FALSE, edge.width = 0.25)

# write full, time calibrated tier 1 tree

finaltreeFile <- 'best_ultrametric_fulltree_ddBD_revision.tre'
write.tree(finalFullTree, finaltreeFile)


# How does phytools::force_ultrametric compare to ape::chronoMPL?
## ape::chronoMPL is much much faster than force_ultrametric, and since we will be applying this to all tier2 trees as well, we will use chronoMPL. Differences are negligable.

maxAge <- max(c(branching.times(finalFullTree), branching.times(finalFullTree2)))
plot(branching.times(finalFullTree), branching.times(finalFullTree2), xlim = c(0, maxAge), ylim = c(0, maxAge))
abline(a = 0, b = 1)
cor.test(branching.times(finalFullTree), branching.times(finalFullTree2))
range(abs(branching.times(finalFullTree) - branching.times(finalFullTree2)))


