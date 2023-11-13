## Time calibration of pseudo-posterior

# Time calibrate full tree with treePL, by breaking down tree into a couple of subclades.

## Use the smoothing values identified for main tree.
## STEPS:
## 	- for each unique tier2 tree, pull a sample from MCMCTREE posterior, and write the treePL config file to folder.
## 	- run all treePL runs
## 	- reassemble

library(ape)
library(castor)
library(MCMCtreeR)

backboneFile <- '~/Dropbox/squamatePhylo/2019/part10/withOutgroups/ddbd_ARS3_June2023/mcmctree_AR/calibS3v2/mcmctree_useData2/iterA/FigTree.tre'
mcmcPostFile <- '~/Dropbox/squamatePhylo/2019/part10/withOutgroups/ddbd_ARS3_June2023/mcmctree_AR/calibS3v2/mcmctree_useData2/iterA/mcmc.txt'

fossilfile <- '~/Dropbox/squamatePhylo/2019/squamate_fossil_calibrations_v6.csv'



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

backboneTree <- readMCMCtree(backboneFile)
mcmcPost <- data.table::fread(mcmcPostFile, data.table = FALSE)

subtreeTips <- list(
    
    gekko = c('Coleonyx_variegatus', 'Crenadactylus_naso', 'Nephrurus_levis', 'Lialis_burtonis', 'Hemidactylus_frenatus', 'Gonatodes_humeralis', 'Thecadactylus_solimoensis'),
    
    scinco = c('Plestiodon_skiltonianus', 'Zonosaurus_laticaudatus', 'Smaug_mossambicus', 'Xantusia_vigilis'),
    
    laterata = c('Bipes_biporus', 'Monopeltis_guentheri', 'Dryadosaura_nordestina', 'Ameiva_ameiva', 'Lacerta_viridis', 'Rhineura_floridana', 'Ptychoglossus_brevifrontalis', 'Blanus_strauchi', 'Cadea_blanoides', 'Trogonophis_wiegmanni'),
    
    iguania = c('Chamaeleo_calyptratus', 'Moloch_horridus', 'Uma_notata', 'Dipsosaurus_dorsalis', 'Anolis_punctatus', 'Diplolaemus_darwinii', 'Polychrus_liogaster', 'Hoplocercus_spinosus', 'Leiocephalus_macropus', 'Phymaturus_palluma', 'Oplurus_cyclurus', 'Basiliscus_basiliscus', 'Crotaphytus_collaris', 'Stenocercus_squarrosus'),
    
    anguids = c('Heloderma_suspectum', 'Ophiodes_striatus', 'Anniella_pulchra', 'Elgaria_kingii', 'Xenosaurus_platyceps', 'Shinisaurus_crocodilurus', 'Varanus_griseus', 'Lanthanotus_borneensis'),
    
    snakes = c('Acrochordus_granulatus', 'Anilius_scytale', 'Liotyphlops_ternetzii', 'Charina_bottae', 'Casarea_dussumieri', 'Contia_tenuis', 'Cylindrophis_ruffus', 'Micrurus_fulvius', 'Gerrhopilus_mirus', 'Fordonia_leucobalia', 'Atractaspis_fallax', 'Namibiana_labialis', 'Loxocemus_bicolor', 'Python_regius', 'Trachyboa_boulengeri', 'Anilios_waitii', 'Rhinophis_saffragamus', 'Agkistrodon_contortrix', 'Xenopeltis_unicolor', 'Anomochilus_leonardi', 'Pareas_nuchalis', 'Kladirostratus_acutus', 'Fimbrios_klossi', 'Xenophidion_schaeferi', 'Xenotyphlops_grandidieri')
    
)


# What were the optimal smoothing parameters inferred for tier1
tier1SmoothParams <- numeric(length(subtreeTips) + 1)
names(tier1SmoothParams) <- c(names(subtreeTips), 'scaffold')
tier1SmoothParams['gekko'] <- 1
tier1SmoothParams['scinco'] <- 1
tier1SmoothParams['laterata'] <- 1
tier1SmoothParams['iguania'] <- 0.1
tier1SmoothParams['anguids'] <- 0.1
tier1SmoothParams['snakes'] <- 0.1
tier1SmoothParams['scaffold'] <- 1

tier1optParams <- numeric(length(subtreeTips) + 1)
names(tier1optParams) <- c(names(subtreeTips), 'scaffold')
tier1optParams['gekko'] <- 554
tier1optParams['scinco'] <- 552
tier1optParams['laterata'] <- 554
tier1optParams['iguania'] <- 555
tier1optParams['anguids'] <- 213
tier1optParams['snakes'] <- 552
tier1optParams['scaffold'] <- 211




# Read in tier 2 trees

setwd('~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/tier2/fulltrees/finaltrees')

tier2trees <- list.files(pattern = '\\.tre$')
tier2trees <- lapply(tier2trees, read.tree)

tier2trees <- lapply(tier2trees, function(x) root(x, 'Sphenodon_punctatus', resolve.root = TRUE))
tier2trees <- lapply(tier2trees, function(x) read.tree(text = write.tree(ladderize(x), file = '')))
tier2trees <- lapply(tier2trees, function(x) drop.tip(x, 'Sphenodon_punctatus'))

class(tier2trees) <- 'multiPhylo'

# do any tier2 trees have the same topology as backbone tree?
table(sapply(tier2trees, function(x) all(backboneTree[[1]]$tip.label %in% x$tip.label)))
# sapply(tier2trees, function(x) all.equal.phylo(backboneTree[[1]], keep.tip(x, backboneTree[[1]]$tip.label), use.edge.length = FALSE))
# none do.

phangorn::RF.dist(tier2trees)

# read in file that shows what percentage of the bootstrap distribution each unique topology represents
percentTopo <- read.csv('~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/tier2/percentTopologyInBootstrapPosterior.csv')[,2]
names(percentTopo) <- paste0('topo', 1:29)
sort(percentTopo)
as.matrix(percentTopo)
percentTopo[percentTopo >= 0.01]

length(unique(sample(names(percentTopo), size = 100, replace = TRUE, prob = percentTopo)))

topoRep <- replicate(n = 1000, sample(names(percentTopo), size = 100, replace = TRUE, prob = percentTopo), simplify = FALSE)
sort(sapply(topoRep, function(x) length(unique(x))))
sort(table(sapply(topoRep, function(x) length(unique(x)))))

set.seed(100)
tier2trees100 <- sample(names(percentTopo), size = 100, replace = TRUE, prob = percentTopo)
length(unique(tier2trees100))

# are there any tip name issues?
setdiff(backboneTree[[1]]$tip.label, tier2trees[[1]]$tip.label) # just outgroups, good!

commonSp <- intersect(tier2trees[[1]]$tip.label, backboneTree[[1]]$tip.label)

prunedTier2 <- lapply(tier2trees, function(x) keep.tip(x, commonSp))
class(prunedTier2) <- 'multiPhylo'

phangorn::RF.dist(prunedTier2)
# all tier2 trees are identical when reduced to taxa in calibrated backbone

# For each internal node in calibrated backbone, check that the same node is present in tier 2 trees
# (ignore nodes that have tips not found at all in tier2 trees (outgroups))
## nodeCheck matrix contains nodes as rows, tier2 trees as columns, and presence/absence of equivalent node

backboneNodes <- (Ntip(backboneTree[[1]]) + 1):Nnode(backboneTree[[1]], internal.only = FALSE)
nodeCheck <- matrix(nrow = length(backboneNodes), ncol = length(prunedTier2))
for (i in 1:length(backboneNodes)) {
	
	backboneClade <- extract.clade(backboneTree[[1]], backboneNodes[i])
	if (all(backboneClade$tip.label %in% commonSp)) {
		nodeCheck[i, ] <- sapply(prunedTier2, function(x) identical(sort(backboneClade$tip.label), sort(extract.clade(x, getMRCA(x, backboneClade$tip.label))$tip.label)))
		message('\t', i, '-- ', sum(sapply(prunedTier2, function(x) identical(sort(backboneClade$tip.label), sort(extract.clade(x, getMRCA(x, backboneClade$tip.label))$tip.label)))))
	} else {
		# for nodes where taxa are not all present in tree (due to outgroups), mark as TRUE, just so that this doesn't interfere
		nodeCheck[i, ] <- TRUE
	}
}

# nodeCheck: each row is a node, each column is a tree
table(apply(nodeCheck, 2, function(x) sum(x, na.rm = TRUE)))
table(apply(nodeCheck, 1, function(x) sum(x, na.rm = TRUE)))
xx <- which(apply(nodeCheck, 1, function(x) sum(x, na.rm = TRUE)) != 29)

# There are 2 nodes that not identical in tier2 (in any tier2 tree)

spSet <- Reduce(union, lapply(prunedTier2, function(x) extract.clade(x, getMRCA(x, extract.clade(backboneTree[[1]], backboneNodes[xx[1]])$tip.label))$tip.label))
treeSet <- lapply(prunedTier2, function(x) keep.tip(x, spSet))
class(treeSet) <- 'multiPhylo'
phangorn::RF.dist(treeSet)

spSet2 <- Reduce(union, lapply(prunedTier2, function(x) extract.clade(x, getMRCA(x, extract.clade(backboneTree[[1]], backboneNodes[xx[2]])$tip.label))$tip.label))
treeSet2 <- lapply(prunedTier2, function(x) keep.tip(x, spSet2))
class(treeSet2) <- 'multiPhylo'
phangorn::RF.dist(treeSet2)


excludeNodes <- backboneNodes[xx]


## We will exclude these two nodes when time calibrating tier2 trees

fossildat <- read.csv(fossilfile, stringsAsFactors=FALSE)

fossildat <- fossildat[which(fossildat[,1] != ''),]
fossildat <- fossildat[which(fossildat$use != ''),]
fossildat <- fossildat[which(fossildat$use == 'yes'),]
fossildat <- fossildat[which(fossildat$tree == 'genomicBackbone'),]


# subset table to only fossils for which spanning taxa are in tree
useFossils <- c()

for (i in 1:nrow(fossildat)) {
	spanningTaxa <- strsplit(fossildat[i, 'spanning.taxa'], split = '\\|')[[1]]
	if (all(spanningTaxa %in% tier2trees[[1]]$tip.label)) {
		useFossils <- c(useFossils, i)
	}
}

# not in tree, but otherwise flagged for use
fossilsNotUsed <- fossildat[setdiff(1:nrow(fossildat), useFossils), ]

fossildat <- fossildat[useFossils, ]

# prep all unique combinations of calibrated tree tip labels
zz <- backboneTree[[1]]$tip.label
calibTipPairs <- expand.grid(list(zz, zz))
calibTipPairs[,1] <- as.character(calibTipPairs[,1])
calibTipPairs[,2] <- as.character(calibTipPairs[,2])
calibTipPairs <- lapply(1:nrow(calibTipPairs), function(i) sort(as.character(calibTipPairs[i, ])))
calibTipPairs <- do.call(rbind, calibTipPairs)
calibTipPairs <- calibTipPairs[which(apply(calibTipPairs, 1, function(x) x[1] == x[2]) == FALSE), ]
calibTipPairs <- unique(calibTipPairs)
calibTipPairs <- calibTipPairs[apply(calibTipPairs, 1, function(x) all(x %in% tier2trees[[1]]$tip.label)), ]


# for each unique tier2 genomic topology, we will sample 100 mcmctree posteriors, x the number of times this topology appears in bootstrap posterior out of 100.

# Create a directory for each, named topoXXrunYY
maindir <- '~/Dropbox/squamatePhylo/2019/treePL_fulltree/resubmission/tier2'

# test
totalTrees <- 0
for (i in 1:length(tier2trees)) {
    # representation in set of trees * number of posterior mcmctree samples we will use
    nReps <- sum(tier2trees100 == paste0('topo', i)) * 1
    
    # totalTrees <- totalTrees + nReps
    
    # sample from mcmctree posterior
    mcmctreeSamples <- sample(1:nrow(mcmcPost), size = nReps, replace = FALSE)
    totalTrees <- totalTrees + length(mcmctreeSamples)
    
    
}
totalTrees

reps <- paste0('iter', 1:10)

for (i in 1:length(tier2trees)) {
	
	# representation in set of trees * number of posterior mcmctree samples we will use
	nReps <- sum(tier2trees100 == paste0('topo', i)) * 1

	if (nReps > 0) {		
	    # match up nodes in mcmcPost to nodes in full tree
	    mcmcCols <- grep('t_n', colnames(mcmcPost), value = TRUE)
	    mcmcCorrespondingNodes <- c()
	    for (j in 1:length(mcmcCols)) {
	        if (all(geiger::tips(backboneTree[[1]], as.numeric(gsub('t_n', '', mcmcCols[j]))) %in% tier2trees[[i]]$tip.label)) {
	            # tier2mrca <- getMRCA(tier2trees[[i]], geiger::tips(backboneTree[[1]], as.numeric(gsub('t_n', '', mcmcCols[j]))))
	            # backboneClade <- extract.clade(backboneTree[[1]], as.numeric(gsub('t_n', '', mcmcCols[j])))
	            # tier2Clade <- keep.tip(extract.clade(tier2trees[[i]], tier2mrca), backboneClade$tip.label)
	            # all.equal.phylo(backboneClade, tier2Clade, use.edge.length = FALSE)
	           
	            mcmcCorrespondingNodes[j] <- getMRCA(tier2trees[[i]], geiger::tips(backboneTree[[1]], as.numeric(gsub('t_n', '', mcmcCols[j]))))
	        } else {
	            mcmcCorrespondingNodes[j] <- NA
	        }
	    }
	    
	    # drop those that don't have an equivalent node in tier2 trees (due to outgroups)
	    mcmcCols <- mcmcCols[!is.na(mcmcCorrespondingNodes)]
	    mcmcCorrespondingNodes <- mcmcCorrespondingNodes[!is.na(mcmcCorrespondingNodes)]
	    
	    # drop the nodes that did not have equivalent clades in the tier2 trees
	    mcmcCorrespondingNodes <- mcmcCorrespondingNodes[!mcmcCols %in% paste0('t_n', excludeNodes)]
	    mcmcCols <- mcmcCols[!mcmcCols %in% paste0('t_n', excludeNodes)]
	    
	    length(mcmcCorrespondingNodes) == length(mcmcCols)
	    
	    # sample from mcmctree posterior
	    mcmctreeSamples <- sample(1:nrow(mcmcPost), size = nReps, replace = FALSE)
	    
	    # allMRCA <- apply(calibTipPairs, 1, function(x) getMRCA(tier2trees[[i]], x))
	    
	    for (j in 1:length(mcmctreeSamples)) {
	        calibdat <- as.data.frame(matrix(nrow = length(mcmcCols), ncol = 6))
	        colnames(calibdat) <- c('MinAge', 'MaxAge', 'fossil', 'fulltreeNode', 'backboneNode', 'spanningTaxa')	
	        
	        for (k in 1:length(mcmcCols)) {
	            
	            calibdat[k, 'backboneNode'] <- as.numeric(gsub('t_n', '', mcmcCols[k]))
	            calibdat[k, 'fulltreeNode'] <- mcmcCorrespondingNodes[k]
	            calibdat[k, 'fossil'] <- 'secondaryCalibration'
	            if (is.na(mcmcCorrespondingNodes[k])) {
	                calibdat[k, 'spanningTaxa'] <- NA
	            } else {	
	                calibdat[k, 'spanningTaxa'] <- paste0(getSpanningTips(tier2trees[[i]], mcmcCorrespondingNodes[k]), collapse = ' ')
	                #calibdat[k, 'spanningTaxa'] <- calibTipPairs[which(allMRCA == mcmcCorrespondingNodes[k])[1], ]
	           }
	            calibdat[k, c('MinAge', 'MaxAge')] <- mcmcPost[mcmctreeSamples[j], mcmcCols[k]] * 100
	            
	        }
	        
	        calibdat$MinAge <- as.numeric(calibdat$MinAge)
	        calibdat$MaxAge <- as.numeric(calibdat$MaxAge)
	        
	        treedir <- paste0(maindir, '/genomTopo', i, '_rep', j)
	        if (!dir.exists(treedir)) dir.create(treedir)
	        
	        # for each subtree, set up treePL run
	        for (k in 1:length(subtreeTips)) {
	            
	            cladeName <- names(subtreeTips)[k]
	            
	            subtreedir <- paste0(treedir, '/', cladeName)
	            if (!dir.exists(subtreedir)) dir.create(subtreedir)
	            setwd(subtreedir)
	            
	            cladeTree <- extract.clade(tier2trees[[i]], getMRCA(tier2trees[[i]], subtreeTips[[k]]))
	            cladeTreeName <- paste0(cladeName, '.tre')
	            write.tree(cladeTree, cladeTreeName)
	            
	            # determine which calibrations will be used
	            inClade <- c()
	            for (l in 1:nrow(calibdat)) {
	                inClade[l] <- all(strsplit(calibdat[l, 'spanningTaxa'], ' ')[[1]] %in% cladeTree$tip.label)
	            }
	            
	            # write config file for treePL

				for (l in 1:length(reps)) {	            
		            treefile <- cladeTreeName
		            nthreads <- 1
		            numsites <- 72152
		            
		            treePLfile <- paste0('treePL_config_tier2_fulltree_genomTopo', i, '_rep', j, '_', reps[l], '.txt')
		            outfileTag <- paste0(cladeName, '_treePL_', reps[l])
		            
		            write('# treePL config file for time calibration of tier2 phylogeny.', file = treePLfile, append = FALSE)
		            write(paste('# written', Sys.time()), file = treePLfile, append = TRUE)
		            
		            write(paste0('\n\ntreefile = ', treefile), file = treePLfile, append = TRUE)
		            write(paste0('nthreads = ', nthreads), file = treePLfile, append = TRUE)
		            write(paste0('numsites = ', numsites), file = treePLfile, append = TRUE)
		            
		            write('\n###############################\n', file = treePLfile, append = TRUE)
		            
		            for (l in 1:nrow(calibdat)) {
		                
		                if (inClade[l]) {
		                    
		                    write(paste0('mrca = calibNode_', calibdat[l, 'fulltreeNode'], ' ', calibdat[l, 'spanningTaxa']), file = treePLfile, append = TRUE)
		                    write(paste0('min = calibNode_', calibdat[l, 'fulltreeNode'], ' ', round(calibdat[l, 'MinAge'], 4)), file = treePLfile, append = TRUE)
		                    write(paste0('max = calibNode_', calibdat[l, 'fulltreeNode'], ' ', round(calibdat[l, 'MaxAge'], 4)), file = treePLfile, append = TRUE)
		                    write('', file = treePLfile, append = TRUE)
		                }
		            }
		            
		            write('###############################', file = treePLfile, append = TRUE)
		            
		            write('\n###############################\n', file = treePLfile, append = TRUE)
		            
		            # optimized params
		            treePLopt <- as.numeric(strsplit(as.character(tier1optParams[cladeName]), '')[[1]])
		            
		            # initial parameter block for verifying that different options for opt/optad/optcvad identify the same smooth value
		            write('thorough', file = treePLfile, append = TRUE)
		            write(paste0('opt = ', treePLopt[1]), file = treePLfile, append = TRUE)
		            write(paste0('optad = ', treePLopt[2]), file = treePLfile, append = TRUE)
		            write(paste0('optcvad = ', treePLopt[3]), file = treePLfile, append = TRUE)
		            
		            # smoothing parameter as determined for tier1
		            write(paste0('\nsmooth = ', tier1SmoothParams[cladeName]), file = treePLfile, append = TRUE)
		            		            
		            write(paste0('\noutfile = ', outfileTag, '.tre'), file = treePLfile, append = TRUE)
				}
		            
	        }
	        
	        # scaffold
	        
	        subtreedir <- paste0(treedir, '/scaffoldTree')
	        if (!dir.exists(subtreedir)) dir.create(subtreedir)
	        setwd(subtreedir)
	        
	        # create folder and files for tree that will act as scaffolding for subtrees. This tree will contain the tips used to define the subtrees, PLUS all tips found outside these subtrees. This way, each subtree can be grafted. 
	        
	        otherTaxa <- setdiff(tier2trees[[i]]$tip.label, unlist(lapply(subtreeTips, function(x) extract.clade(tier2trees[[i]], getMRCA(tier2trees[[i]], x))$tip.label)))
	        cladeTree <- keep.tip(tier2trees[[i]], c(otherTaxa, unlist(subtreeTips)))
	        
	        any(otherTaxa %in% backboneTree[[1]]$tip.label)
	        
	        # We will not include secondary calibrations from within the subtrees, since those would be arbitrarily included based on taxon choice
	        
	        # get taxa that are in scaffold tree and that define crown nodes
	        scaffoldNodeTaxa <- lapply(subtreeTips, function(x) {
	            
	            allpairs <- expand.grid(list(x, x))
	            allpairs[,1] <- as.character(allpairs[,1])
	            allpairs[,2] <- as.character(allpairs[,2])
	            allpairs <- apply(allpairs, 1, function(y) sort(y), simplify = FALSE)
	            # allpairs <- lapply(1:nrow(allpairs), function(y) as.character(sort(allpairs[y, ])))
	            allpairs <- allpairs[which(sapply(allpairs, function(y) y[1] == y[2]) == FALSE)]
	            allpairs <- unique(allpairs)
	            
	            allpairs[[which(sapply(allpairs, function(y) getMRCA(tier2trees[[i]], y)) == getMRCA(tier2trees[[i]], x))[1]]]
	        })
	        
	        scaffoldNodeTaxa <- unlist(c(scaffoldNodeTaxa, otherTaxa))
	        
	        # determine what node numbers this tree would contain in the full tree
	        allpairs <- expand.grid(list(scaffoldNodeTaxa, scaffoldNodeTaxa))
	        allpairs[,1] <- as.character(allpairs[,1])
	        allpairs[,2] <- as.character(allpairs[,2])
	        allpairs <- apply(allpairs, 1, function(y) sort(y), simplify = FALSE)
	        # allpairs <- lapply(1:nrow(allpairs), function(x) as.character(sort(allpairs[x, ])))
	        allpairs <- allpairs[which(sapply(allpairs, function(x) x[1] == x[2]) == FALSE)]
	        allpairs <- unique(allpairs)
	        scaffoldNodes <- sapply(allpairs, function(x) getMRCA(tier2trees[[i]], x))
	        
	        inClade <- calibdat$fulltreeNode %in% scaffoldNodes
	        sum(inClade)
	        
	        # plot(cladeTree, cex = 0.4)
	        # calibratedScaffoldNodes <- unique(sapply(allpairs[which((scaffoldNodes %in% calibdat$fulltreeNode) == TRUE)], function(x) getMRCA(cladeTree, x)))
	        # nodelabels(rep('', length(calibratedScaffoldNodes)), node = calibratedScaffoldNodes, frame = 'circle', cex = 0.3, bg = 'orange')
	        
	        cladeTreeName <- 'scaffold.tre'
	        write.tree(cladeTree, cladeTreeName)
	        
	        # write config file for treePL
	        
	        treefile <- cladeTreeName
	        nthreads <- 1
	        numsites <- 72152

			for (l in 1:length(reps)) {	        
		        treePLfile <- paste0('treePL_config_tier2_scaffold_genomTopo', i, '_rep', j, '_', reps[l],  '.txt')
		        outfileTag <- paste0('scaffold_treePL_', reps[l])
		        
		        write('# treePL config file for time calibration of tier2 phylogeny.', file = treePLfile, append = FALSE)
		        write(paste('# written', Sys.time()), file = treePLfile, append = TRUE)
		        
		        write(paste0('\n\ntreefile = ', treefile), file = treePLfile, append = TRUE)
		        write(paste0('nthreads = ', nthreads), file = treePLfile, append = TRUE)
		        write(paste0('numsites = ', numsites), file = treePLfile, append = TRUE)
		        
		        write('\n###############################\n', file = treePLfile, append = TRUE)
		        
		        for (l in 1:nrow(calibdat)) {
		            
		            if (inClade[l]) {
		                
		                spanningpair <- allpairs[[which(scaffoldNodes == calibdat[l, 'fulltreeNode'])[1]]]
		                spanningpair <- paste(spanningpair, collapse = ' ')
		                
		                write(paste0('mrca = calibNode_', calibdat[l, 'fulltreeNode'], ' ', spanningpair), file = treePLfile, append = TRUE)
		                write(paste0('min = calibNode_', calibdat[l, 'fulltreeNode'], ' ', round(calibdat[l, 'MinAge'], 4)), file = treePLfile, append = TRUE)
		                write(paste0('max = calibNode_', calibdat[l, 'fulltreeNode'], ' ', round(calibdat[l, 'MaxAge'], 4)), file = treePLfile, append = TRUE)
		                write('', file = treePLfile, append = TRUE)
		                
		            }
		        }
		        
		        write('###############################', file = treePLfile, append = TRUE)
		        
		        write('\n###############################\n', file = treePLfile, append = TRUE)
		        
		        # optimized params
		        treePLopt <- as.numeric(strsplit(as.character(tier1optParams['scaffold']), '')[[1]])
		        
		        # initial parameter block for verifying that different options for opt/optad/optcvad identify the same smooth value
		        write('thorough', file = treePLfile, append = TRUE)
		        write(paste0('opt = ', treePLopt[1]), file = treePLfile, append = TRUE)
		        write(paste0('optad = ', treePLopt[2]), file = treePLfile, append = TRUE)
		        write(paste0('optcvad = ', treePLopt[3]), file = treePLfile, append = TRUE)
		        
		        # smoothing parameter as determined for tier1
		        write(paste0('\nsmooth = ', tier1SmoothParams['scaffold']), file = treePLfile, append = TRUE)
		        
		        write(paste0('\noutfile = ', outfileTag, '.tre'), file = treePLfile, append = TRUE)
			}	        
	    }
	}
}

# prepare execute file
setwd(maindir)
dirnames <- list.files(pattern = 'genomTopo')
dirnames <- dirnames[dir.exists(dirnames)]

allRuns <- c()
for (i in 1:length(dirnames)) {

    setwd(dirnames[i])
    subdirs <- list.files()    
    subdirs <- subdirs[dir.exists(subdirs)]
    
    configFiles <- list.files(pattern = 'treePL_config', recursive = TRUE)
 
    # prep path and call
    outputFileNames <- paste0('output_iter', gsub('(.+)(iter)(\\d\\d?)(\\.txt$)', '\\3', configFiles), '.txt')
    cmd <- paste0('cd ', dirnames[i], '/', dirname(configFiles), '; treePL ', basename(configFiles), ' > ', outputFileNames, ';')
 
    allRuns <- c(allRuns, cmd)

    setwd(maindir)
}

write(allRuns, file = 'runTier2treePL.txt')

## run with parallel < runTier2treePL.txt --bar --joblog tier2runs.txt
### the joblog will list which jobs have completed, and their runtimes.



## Move through directories and identify the best runs

setwd(maindir)
dirnames <- list.files(pattern = 'genomTopo')
dirnames <- dirnames[dir.exists(dirnames)]

problems <- c()
for (i in 1:length(dirnames)) {
	
	setwd(dirnames[i])
    message('\t', dirnames[i])
	
	subfolders <- list.files()
	subfolders <- subfolders[dir.exists(subfolders)]
	
	for (j in 1:length(subfolders)) {
		setwd(subfolders[j])

		ff <- list.files(pattern = 'output')
	    ll <-  numeric(length(ff))
	    names(ll) <- gsub('(output_)(iter\\d+)\\.txt', '\\2', ff)
	    for (k in 1:length(ff)) {
	        out <- scan(ff[k], what = 'character', sep = '\n', quiet = TRUE)
	        if (grepl('(after opt calc:\\s)(\\d+(\\.\\d+)?)', tail(out, 1))) {
	        	score <- as.numeric(gsub('(after opt calc:\\s)(\\d+(\\.\\d+)?)', '\\2', tail(out, 1)))
	        } else {
	        	score <- NA
	        }
	        
	        ll[k] <- score
	    }
	    # plot(sort(ll))
	    if (!all(is.na(ll))) {
	    	if (anyNA(ll)) stop()
		    treeFile <- grep('iter\\d\\d?\\.tre$', list.files(), value = TRUE)
		    treeFile <- grep(paste0(names(ll)[which.min(ll)], '\\.'), treeFile, value = TRUE)
			file.copy(from = treeFile, to = gsub('iter\\d\\d?', 'best', treeFile))	
		} else {
			problems <- c(problems, paste0(dirnames[i], '_', subfolders[j]))	
		}
		setwd(maindir)
		setwd(dirnames[i])
	}
	setwd(maindir)
}


## Read in tier2 subtrees and assemble into full trees

for (i in 1:length(dirnames)) {

    setwd(dirnames[i])
    
    message('\t', dirnames[i])

    finalsubtrees <- vector('list', length(subtreeTips))
    names(finalsubtrees) <- names(subtreeTips)
    for (j in 1:length(subtreeTips)) {
        finalsubtrees[[j]] <- try(read.tree(list.files(path = names(subtreeTips)[j], pattern = 'treePL_best\\.tre$', full.names = TRUE)), silent = TRUE)
    }
    
    if (any(sapply(finalsubtrees, class) == 'NULL')) {
        setwd(maindir)
        next
    }
    
    
    sapply(finalsubtrees, function(x) is.binary(x))
    sapply(finalsubtrees, function(x) is.ultrametric(x, option = 2))
    sapply(finalsubtrees, function(x) Ntip(x))

    # lapply(finalsubtrees, function(x) keep.tip(x, sample(intersect(x$tip.label, finalScaffold$tip.label), 3)))

    finalScaffold <- read.tree('scaffoldTree/scaffold_treePL_best.tre')
    is.binary(finalScaffold)
    is.ultrametric(finalScaffold, option = 2)

    # add tag to more easily be able to remove later
    finalScaffold$tip.label <- paste0(finalScaffold$tip.label, '_scaffold')

    # assemble
    finalFullTree <- finalScaffold

    # pdf('finalScaffoldAnnotated.pdf', width = 20, height = 20)
    # plot.phylo(finalScaffold, cex = 0.7, edge.width = 0.6)
    # nodelabels(round(branching.times(finalScaffold), 2), node = as.numeric(names(branching.times(finalScaffold))), cex = 0.7)
    # dev.off()

    # subtreeTips lists the subset of species that represent each subtree. We will swap out that subset with the full, time-calibrated tree that was calibrated separately. Crown ages should match.
    for (j in 1:length(subtreeTips)) {

        finalScaffoldCopy <- finalScaffold
        finalScaffoldCopy$tip.label <- gsub('_scaffold$', '', finalScaffoldCopy$tip.label)
        branching.times(finalScaffoldCopy)[as.character(getMRCA(finalScaffoldCopy, subtreeTips[[j]]))]
        max(branching.times(finalsubtrees[[j]]))
        # they match. good.


        # replace finalFullTree with one where the subtree has been added in
        finalFullTreeCopy <- finalFullTree
        finalFullTreeCopy$tip.label <- gsub('_scaffold$', '', finalFullTreeCopy$tip.label)
        bindNode <- getMRCA(finalFullTreeCopy, subtreeTips[[j]])
        finalFullTree <- bind.tree(finalFullTree, finalsubtrees[[j]], where = bindNode, position = 0)
        toDrop <- paste0(finalsubtrees[[j]]$tip.label, '_scaffold')[which(paste0(finalsubtrees[[j]]$tip.label, '_scaffold') %in% finalFullTree$tip.label)]
        finalFullTree <- drop.tip(finalFullTree, toDrop)

        # branching time should again match.
        branching.times(finalFullTree)[as.character(getMRCA(finalFullTree, subtreeTips[[j]]))]

        # should be ultrametric still!
        is.ultrametric(finalFullTree, option = 2)
        is.binary(finalFullTree)
    }

    length(unique(finalFullTree$tip.label))
    Ntip(finalFullTree)

    grep('_scaffold$', finalFullTree$tip.label, value = TRUE)
    finalFullTree$tip.label <- gsub('_scaffold$', '', finalFullTree$tip.label)
    grep('_scaffold$', finalFullTree$tip.label, value = TRUE)

    finalFullTree <- ape::chronoMPL(finalFullTree)

    # final tree is not ultrametric, but only because of numerical precision. We will fix this with phytools::force.ultrametric.
    is.ultrametric(finalFullTree)
    is.ultrametric(finalFullTree, option = 2)
    is.binary(finalFullTree)

    # all.equal.phylo(finalFullTree, tier2trees[[as.numeric(gsub('genomTopo(\\d\\d?)_rep\\d\\d?\\d?', '\\1', dirnames[i]))]], use.edge.length = FALSE)

    finaltreesDir <- '../tier2finaltrees'
    if (!dir.exists(finaltreesDir)) dir.create(finaltreesDir)

    write.tree(finalFullTree, paste0(finaltreesDir, '/tier2_ultrametric_', dirnames[i], '.tre'))

    setwd(maindir)
}



# Comparisons to tier1 tree
tier1 <- read.tree('~/Dropbox/Oz_Crown_Ages/final-trees/best_ultrametric_fulltree_ddBD_revision.tre')

tier2 <- vector('list', length = length(list.files('tier2finaltrees', pattern = '\\.tre$')))
for (i in 1:length(list.files('tier2finaltrees', pattern = '\\.tre$'))) {
    tier2[[i]] <- read.tree(list.files('tier2finaltrees', pattern = '\\.tre$', full.names = TRUE)[i])
} 
class(tier2) <- 'multiPhylo'

table(sapply(tier2, is.binary))
table(sapply(tier2, is.ultrametric))

library(phangorn)

# RF dists
rf <- RF.dist(unroot(tier1), unroot(tier2))
table(rf)

names(tier2) <- list.files('tier2finaltrees', pattern = '\\.tre$')

write.tree(tier2, file = 'tier2_100.trees', tree.names = names(tier2))


# do tier2 trees match genomic constraint?
tier1Constraint <- read.tree('~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/tier1/families/fullGenomBackbone_tier1.tre')
tier2 <- read.tree('~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/tier2/tier2_100.trees')
tier2Constraints <- read.tree('~/Dropbox/Oz_Crown_Ages/phylogenetic_inference/tier2/constraintTrees.trees')

all(tier1Constraint$tip.label %in% tier2Constraints[[1]]$tip.label)
tier2Constraints <- lapply(tier2Constraints, function(x) keep.tip(x, intersect(tier1Constraint$tip.label, x$tip.label)))

all(tier2Constraints[[1]]$tip.label %in% tier2[[1]]$tip.label)
setdiff(tier2Constraints[[1]]$tip.label,tier2[[1]]$tip.label) # this is expected: Gehyra_Cysp vs CYsp

tier2Constraints <- lapply(tier2Constraints, function(x) keep.tip(x, intersect(tier2[[1]]$tip.label, x$tip.label)))

problem <- c()
for (i in 1:length(tier2)) {
    # we expect a tier2 tree to match the topology of one of the tier2constraint trees
    ct <- sum(sapply(tier2Constraints, function(x) all.equal.phylo(keep.tip(tier2[[1]], x$tip.label), x, use.edge.length = FALSE)))
    if (ct != 1) problem <- c(problem, i)
}

length(problem)

# For the two nodes that had to be ignored, how does the distribution of divergence times compare to the tier1 divergence time?

spSet1 <- c('Farancia_abacura', 'Dendrelaphis_calligaster')
spSet2 <- c('Stenocercus_squarrosus', 'Liolaemus_valdesianus')

divTimes1 <- numeric(length(tier2))
divTimes2 <- numeric(length(tier2))
for (i in 1:length(tier2)) {
	divTimes1[i] <- branching.times(tier2[[i]])[as.character(getMRCA(tier2[[i]], spSet1))]
	divTimes2[i] <- branching.times(tier2[[i]])[as.character(getMRCA(tier2[[i]], spSet2))]	
}

plot(density(divTimes1))
abline(v=branching.times(tier1)[as.character(getMRCA(tier1, spSet1))])

plot(density(divTimes2))
abline(v=branching.times(tier1)[as.character(getMRCA(tier1, spSet2))])

