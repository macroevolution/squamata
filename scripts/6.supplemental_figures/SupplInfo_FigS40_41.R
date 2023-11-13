# Figures S40 and S41: Comparisons of several timetrees (Tonini 2016, Irisarri 2017, Burbrink 2020, this study)

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

library(ape)
library(MCMCtreeR)

alldat <- read.csv('./data/alldat.csv')


fulltree <- read.tree('./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre')

mcmcPostFile <- './data/2.time_calibration_imputation/MCMCTreeAnalyses/mcmctree_useData2/iterA/mcmc.txt'
mcmcPost <- data.table::fread(mcmcPostFile, data.table = FALSE)

tr1 <- readMCMCtree('./data/2.time_calibration_imputation/MCMCTreeAnalyses/mcmctree_useData2/iterA/FigTree.tre')


timeTree <- tr1[[1]]
timeTree <- ladderize(timeTree)
timeTree <- read.tree(text = write.tree(timeTree, file = ''))

nodeTimes <- cbind.data.frame(oldNodes = rownames(tr1[[2]]), newNodes = NA, tr1[[2]])
nodeTimes$newNodes <- sapply(rownames(tr1[[2]]), function(x) getMRCA(timeTree, geiger::tips(tr1[[1]], x)))

# for each node in tree, calculate the median divergence time from the posterior, and the 95% HPD
for (i in 1:nrow(nodeTimes)) {
	xx <- mcmcPost[, paste0('t_n', nodeTimes[i, 'oldNodes'])]
	nodeTimes[i, 'median'] <- median(xx)
	nodeTimes[i, 'hpd5'] <- HDInterval::hdi(xx, credMass = 0.95)['lower']
	nodeTimes[i, 'hpd95'] <- HDInterval::hdi(xx, credMass = 0.95)['upper']
}

# rescale timeTree with median divergence times
newBL <- numeric(length(timeTree$edge.length))
origBL <- numeric(length(timeTree$edge.length))
for (i in 1:nrow(timeTree$edge)) {
	xx <- nodeTimes[nodeTimes$newNodes %in% timeTree$edge[i, ],]
	if (nrow(xx) == 2) {
		origBL[i] <- abs(diff(xx[, 'mean']))
		newBL[i] <- abs(diff(xx[, 'median']))
	} else if (nrow(xx) == 1) {
		origBL[i] <- xx[, 'mean']
		newBL[i] <- xx[, 'median']	
	} else stop()
}

timeTree$edge.length <- newBL

timeTree <- geiger:::rescale.phylo(timeTree, model = 'depth', depth = max(branching.times(timeTree)) * 100)

nodeTimes[, 3:8] <- nodeTimes[, 3:8] * 100



################################
## Comparison of time trees

# Read in trees to compare
tonini <- read.tree('./phylogenetic_inference/treesForComparison/Tonini2016.tre')
burbrink <- read.tree('./phylogenetic_inference/treesForComparison/Burbrink2020.tre')
irisarri <- read.tree('./phylogenetic_inference/treesForComparison/Irisarri2017_CATGTR-LN-BD-SB_100jacks.chronogram_mean_compCrI.tre')

# convert tip labels to species names in Irisarri tree, and associate with reasonable tips from our tree
require(pdftools)

si <- pdf_text('./phylogenetic_inference/treesForComparison/Irisarri_et_al_Supplement_combined.R2.pdf')[[20]]
si <- strsplit(si, '\n')[[1]][6:114]

mat <- matrix(nrow = length(si), ncol = 4)
mat <- as.data.frame(mat)
colnames(mat) <- c('taxonID', 'species', 'dnaType', 'accession')
for (i in 1:length(si)) {
	xx <- strsplit(si[i], '\\s\\s\\s')[[1]]
	xx <- xx[xx != '']
	xx <- gsub('^\\s+|\\s+$', '', xx)
	mat[i,] <- xx
}

table(irisarri$tip.label %in% mat[,1])

keepTips <- c(geiger::tips(irisarri, getMRCA(irisarri, c('Sphenodon_', 'Anolis_car'))), 'Gallus_gal', 'Alligator_', 'Homo_sapie')

irisarri <- keep.tip(irisarri, keepTips)

all(irisarri$tip.label %in% mat[,1])
irisarri$tip.label <- sapply(irisarri$tip.label, function(x) mat[,2][mat[,1] == x])
irisarri$tip.label <- gsub('\\s+', '_', irisarri$tip.label)
irisarri <- ladderize(irisarri)
names(irisarri$tip.label) <- NULL

setdiff(gsub('(.+)_(.+)', '\\1', irisarri$tip.label), gsub('(.+)_(.+)', '\\1', timeTree$tip.label))


# combine
treeList <- list(timeTree = timeTree, burbrink = burbrink, irisarri = irisarri, tonini = tonini)
treeList <- lapply(treeList, ladderize)

# drop outgroups
outgroups <- c('Sphenodon_punctatus', 'Homo_sapiens', 'Gallus_gallus', 'Alligator_mississippiensis')
for (i in 1:length(treeList)) {
    if (any(outgroups %in% treeList[[i]]$tip.label)) {
        treeList[[i]] <- drop.tip(treeList[[i]], intersect(outgroups, treeList[[i]]$tip.label)) 
    }
}


treeNames <- c('this study', 'Burbrink et al. 2020', 'Irisarri et al. 2017', 'Tonini et al. 2016')


# clades for comparison
# Lepidosauria (can't because we do not have sphenodon)
# Squamata
# Gekkota
# Scincoidea
# Lacertoidea
# Serpentes
# Anguimorpha
# Iguania

taxonList <- list(
	squamata = list(
		ourTree = c('Coleonyx_variegatus', 'Thamnophis_butleri'),
		burbrink = c('Coleonyx_variegatus', 'Thamnophis_marcianus'),
		irisarri = c('Eublepharis_macularius', 'Thamnophis_elegans'),
		tonini = c('Coleonyx_variegatus', 'Thamnophis_butleri')
	),
	gekkota = list(
		ourTree = c('Coleonyx_variegatus', 'Delma_tincta'),
		burbrink = c('Coleonyx_variegatus', 'Delma_borea'),
		irisarri = c('Eublepharis_macularius', 'Tarentola_mauritanica'),
		tonini = c('Amalosia_obscura', 'Aeluroscalabotes_felinus')
	),
	scincoidea = list(
		ourTree = c('Cordylus_rivae', 'Saproscincus_lewisi'),
		burbrink = c('Cordylus_sp', 'Lepidophyma_flavimaculatum'),
		irisarri = c('Saproscincus_basiliscus', 'Scincella_lateralis'),
		tonini = c('Acontias_aurantiacus', 'Chamaesaura_aenea')
	),
	lacertoidea = list(
		ourTree = c('Tupinambis_teguixin', 'Amphisbaena_brasiliana'),
		burbrink = c('Tupinambis_teguixin', 'Amphisbaena_fuliginosa'),
		irisarri = c('Tupinambis_teguixin', 'Podarcis_sp.'),
		tonini = c('Ameiva_aggerecusans', 'Acanthodactylus_aureus')
	),
	anguimorpha = list(
		ourTree = c('Shinisaurus_crocodilurus', 'Abronia_taeniata'),
		burbrink = c('Shinisaurus_crocodilurus', 'Anniella_pulchra'),
		irisarri = c('Elgaria_multicarinata', NA),
		tonini = c('Varanus_acanthurus', 'Shinisaurus_crocodilurus')		
	),
	iguania = list(
		ourTree = c('Chamaeleo_calyptratus', 'Iguana_iguana'),
		burbrink = c('Chamaeleo_calyptratus', 'Phrynosoma_cornutum'),
		irisarri = c('Chamaeleo_chamaeleon','Iguana_iguana'),
		tonini = c('Archaius_tigris', 'Callisaurus_draconoides')			
	),
	serpentes = list(
		ourTree = c('Anilios_waitii', 'Hydrophis_platurus'),
		burbrink = c('Indotyphlops_pushpakumara', 'Lampropeltis_getula'),
		irisarri = c('Boa_constrictor', 'Pantherophis_guttatus'),
		tonini = c('Epictia_columbi', 'Liotyphlops_albirostris')	
	)	
)

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

tree <- tonini
getSpanningTips(tree, getMRCA(tree, intersect(tree$tip.label, c('Cordylus_rivae', 'Cordylus_sp', 'Zonosaurus_ornatus', 'Zonosaurus_laticaudatus', 'Lepidophyma_flavimaculatum', 'Saproscincus_lewisi', 'Xantusia_henshawi', 'Xantusia_vigilis', 'Tiliqua_scincoides'))))


all(unlist(lapply(taxonList, function(x) x$ourTree)) %in% timeTree$tip.label)
all(unlist(lapply(taxonList, function(x) x$burbrink)) %in% burbrink$tip.label)
all(unlist(lapply(taxonList, function(x) x$irisarri)) %in% irisarri$tip.label)
all(unlist(lapply(taxonList, function(x) x$tonini)) %in% tonini$tip.label)


includeNodes <- list()

for (i in 1:length(treeList)) {
	
	includeNodes[[i]] <- sapply(taxonList, function(x) {
		if (!anyNA(unlist(x[[i]]))) {
			getMRCA(treeList[[i]], x[[i]])
		} else {
			NA
		}
	})
}


rts <- sapply(treeList, function(x) max(branching.times(x)))
maxrts <- max(rts)
lim <- cbind(rts - maxrts, rts)


pdf('~/Downloads/figS40.pdf', width = 9, height = 15)

par(mfrow=c(length(treeList), 1), mar = rep(0.5, 4), oma = rep(2, 4)) 

nodeCoords <- replicate(length(treeList), matrix(nrow = length(includeNodes[[1]]), ncol = 2), simplify = FALSE)

for (i in 1:length(treeList)) {
	
	plot.phylo(treeList[[i]], x.lim = lim[i,], edge.width = 0.5, show.tip.label = FALSE)

	
	mtext(treeNames[i], side = 2)
	
	# get plot coordinates of nodes of interest
	nodeCoords[[i]][, 1] <- sapply(includeNodes[[i]], function(x) get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[x])
	nodeCoords[[i]][, 2] <- sapply(includeNodes[[i]], function(x) get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[x])
	
	# convert those coordinates to multiplot coordinates
	nodeCoords[[i]][, 1] <- grconvertX(nodeCoords[[i]][, 1], from = 'user', to = 'ndc')
	nodeCoords[[i]][, 2] <- grconvertY(nodeCoords[[i]][, 2], from = 'user', to = 'ndc')
	
	if (i == 1) axisPhylo(3)
	
}

segments(x0 = lim[nrow(lim), 2] - c(50, 100, 150, 200), x1 = lim[nrow(lim), 2] - c(50, 100, 150, 200), y0 = grconvertY(0.02, from = 'ndc', to = 'user'), y1 = grconvertY(0.97, from = 'ndc', to = 'user'), lty = 3, xpd = NA)

phyloch::axisGeo(GTS = strat2012, unit = c("period"), cex = 0.8)


# convert those multiplot coordinates back to current plot coordinates

newCoords <- nodeCoords
for (i in 1:length(nodeCoords)) {
	newCoords[[i]][, 1] <- grconvertX(nodeCoords[[i]][, 1], from = 'ndc')
	newCoords[[i]][, 2] <- grconvertY(nodeCoords[[i]][, 2], from = 'ndc')	
}


# Plot segments that link equivalent nodes
for (i in 1:nrow(newCoords[[1]])) {
	tmp <- do.call(rbind, sapply(newCoords, function(x) x[i,], simplify = FALSE))
	tmp <- tmp[complete.cases(tmp),]
	for (i in 1:(nrow(tmp) - 1)) {
		segments(tmp[i, 1], tmp[i, 2], tmp[i+1, 1], tmp[i+1, 2], lwd = 1, col = 'red', xpd = NA)
	}
}

# plot node labels
points(do.call(rbind, newCoords), pch = 21, bg = 'orange', cex = 2.5, xpd = NA)
text(do.call(rbind, newCoords), labels = rep(1:nrow(newCoords[[1]]), length(treeList)), font = 2, xpd = NA)
	
dev.off()






######################
# Correlations between divergence times

## reduce trees to genus level

burbrink$node.label <- NULL


# return genus level tree
reduceToGenus <- function(tree) {
	
	uniqueGenera <- unique(gsub('(.+)_(.+)', '\\1', tree$tip.label))
	spGroups <- lapply(uniqueGenera, function(x) grep(paste0('^', x, '_'), tree$tip.label, value = TRUE))
	table(lengths(spGroups))
	keepTips <- vector('list', length(uniqueGenera))
	notMono <- c()
	for (i in 1:length(spGroups)) {
		if (length(spGroups[[i]]) > 1) {
			if (Ntip(extract.clade(tree, getMRCA(tree, spGroups[[i]]))) == length(spGroups[[i]])) {
				keepTips[[i]] <- spGroups[[i]][1]
			} else {
				notMono <- c(notMono, i)
				keepTips[[i]] <- spGroups[[i]][1]
			}
		} else {
			keepTips[[i]] <- spGroups[[i]]
		}
	}
	gentree <- keep.tip(tree, unlist(keepTips))
	gentree$tip.label <- gsub('(.+)_(.+)', '\\1', gentree$tip.label)
	list(gentree, uniqueGenera[notMono])
}

reduceToFamily <- function(tree) {
	
	fams <- sapply(tree$tip.label, function(x) alldat[alldat$genus == gsub('(.+)_(.+)', '\\1', x), 'family'][1])
	
	extraGenera <- c(Macropisthodon = 'Colubridae', Seminatrix = 'Colubridae', Pelamis = 'Elapidae', Leiopython = 'Pythonidae', Apodora = 'Pythonidae', Rachodactylus = 'Diplodactylidae', Ptychozoon = 'Gekkonidae', Dryocalamus = 'Colubridae', Colopus = 'Gekkonidae', Androngo = 'Scincidae', Afroablepharus = 'Scincidae', Bassiana = 'Scincidae', Niveoscincus = 'Scincidae', Lepidothyris = 'Scincidae', Hemisphaeriodon = 'Scincidae', Psilophthalmus = 'Gymnophthalmidae', Coloptychon = 'Anguidae', Lemuriatyphlops = 'Typhlopidae', Brachyophidium = 'Uropeltidae', Pseudotyphlops = 'Typhlopidae', Rhagerhis = 'Lamprophiidae', Balanophis = 'Colubridae', Sibynomorphus = 'Colubridae', Chilomeniscus = 'Colubridae', Chionactis = 'Colubridae', Pseustes = 'Colubridae', Cyclophiops = 'Colubridae', Rhinechis = 'Colubridae', Orthriophis = 'Colubridae')
	
	fams[gsub('(.+)_(.+)', '\\1', names(fams)) %in% names(extraGenera)] <- extraGenera[gsub('(.+)_(.+)', '\\1', names(fams)[gsub('(.+)_(.+)', '\\1', names(fams)) %in% names(extraGenera)])]
	
	if (anyNA(fams)) stop()	
	
	uniqueFams <- unique(fams)
	
	keepTips <- sapply(uniqueFams, function(x) names(fams[fams == x][1]))
	
	famTree <- keep.tip(tree, keepTips)
	famTree$tip.label <- fams[famTree$tip.label]
	names(famTree$tip.label) <- NULL
	famTree
	
}

# burbrink
burbrink_genus <- reduceToGenus(burbrink)


# our tree
fulltree_genus <- reduceToGenus(fulltree)

# Tonini
tonini_genus <- reduceToGenus(tonini)

# Irisarri
irisarri_genus <- reduceToGenus(irisarri)

# all non-monophyletic
Reduce(union, list(burbrink_genus[[2]], fulltree_genus[[2]], tonini_genus[[2]], irisarri_genus[[2]]))

# our tree vs Burbrink
xx <- phytools::matchNodes(fulltree_genus[[1]], burbrink_genus[[1]], method = 'descendants')
xx <- xx[complete.cases(xx),]

plot(branching.times(fulltree_genus[[1]])[as.character(xx[,1])], branching.times(burbrink_genus[[1]])[as.character(xx[,2])])


# our tree vs Tonini
fulltree_genus2 <- drop.tip(fulltree_genus[[1]], fulltree_genus[[2]])
tonini_genus2 <- drop.tip(tonini_genus[[1]], tonini_genus[[2]])
xx <- phytools::matchNodes(fulltree_genus2, tonini_genus2, method = 'descendants')
xx <- xx[complete.cases(xx),]

plot(branching.times(fulltree_genus2)[as.character(xx[,1])], branching.times(tonini_genus2)[as.character(xx[,2])])

# family level

for (i in 1:length(treeList)) {
	treeList[[i]]$node.label <- NULL
}

treeList[[1]] <- fulltree

# burbrink
burbrink_fam <- reduceToFamily(treeList[['burbrink']])


# our tree
fulltree_fam <- reduceToFamily(treeList[['timeTree']])

# Tonini
tonini_fam <- reduceToFamily(treeList[['tonini']])

# Irisarri
irisarri_fam <- reduceToFamily(treeList[['irisarri']])


matchNodes <- function(tr1, tr2) {
	
	commonsp <- intersect(tr1$tip.label, tr2$tip.label)
	subtree1 <- keep.tip(tr1, commonsp)
	subtree2 <- keep.tip(tr2, commonsp)
	
	xx <- phytools::matchNodes(subtree1, subtree2, method = 'descendants')
	xx <- xx[complete.cases(xx),]
	
	xx[,1] <- sapply(xx[,1], function(x) getMRCA(tr1, geiger::tips(subtree1, x)))
	xx[,2] <- sapply(xx[,2], function(x) getMRCA(tr2, geiger::tips(subtree2, x)))
	
	return(xx)
}


pdf('~/Downloads/figS41.pdf', width = 12, height = 4)
par(mfrow = c(1,3))

# our tree vs Burbrink
xx <- matchNodes(fulltree_fam, burbrink_fam)

plot.new()
plot.window(xlim = range(c(branching.times(fulltree_fam), branching.times(burbrink_fam))), ylim = range(c(branching.times(fulltree_fam), branching.times(burbrink_fam))))

axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")
mtext('this study', side = 1, line = 2.5, cex = 0.75)
mtext('burbrink et al. 2020', side = 2, line = 2.5, cex = 0.75)

abline(a = 0, b = 1, lty = 2)

points(branching.times(fulltree_fam)[as.character(xx[,1])], branching.times(burbrink_fam)[as.character(xx[,2])], pch = 21, bg = gray(0.85), cex = 1.5)

correlation <- cor.test(branching.times(fulltree_fam)[as.character(xx[,1])], branching.times(burbrink_fam)[as.character(xx[,2])])

# mtext(bquote(italic("r") == .(round(correlation$estimate, 2))), side = 3, line = 1.2)
text(grconvertX(0.90, from = 'npc', to = 'user'), grconvertY(0.05, from = 'npc', to = 'user'), bquote(italic("r") == .(round(correlation$estimate, 2))), xpd = NA, cex = 1.5)

mtext('Burbrink et al. 2020', side = 3, font = 2)


# our tree vs Tonini
xx <- matchNodes(fulltree_fam, tonini_fam)

plot.new()
plot.window(xlim = range(c(branching.times(fulltree_fam), branching.times(tonini_fam))), ylim = range(c(branching.times(fulltree_fam), branching.times(tonini_fam))))

axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")
mtext('this study', side = 1, line = 2.5, cex = 0.75)
mtext('Tonini et al. 2016', side = 2, line = 2.5, cex = 0.75)

abline(a = 0, b = 1, lty = 2)

points(branching.times(fulltree_fam)[as.character(xx[,1])], branching.times(tonini_fam)[as.character(xx[,2])], pch = 21, bg = gray(0.85), cex = 1.5)

correlation <- cor.test(branching.times(fulltree_fam)[as.character(xx[,1])], branching.times(tonini_fam)[as.character(xx[,2])])

# mtext(bquote(italic("r") == .(round(correlation$estimate, 2))), side = 3, line = 1.2)
text(grconvertX(0.90, from = 'npc', to = 'user'), grconvertY(0.05, from = 'npc', to = 'user'), bquote(italic("r") == .(round(correlation$estimate, 2))), xpd = NA, cex = 1.5)

mtext('Tonini et al. 2016', side = 3, font = 2)


# out tree vs Irisarri
xx <- matchNodes(fulltree_fam, irisarri_fam)


plot.new()
plot.window(xlim = range(c(branching.times(fulltree_fam), branching.times(irisarri_fam))), ylim = range(c(branching.times(fulltree_fam), branching.times(irisarri_fam))))

axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1)
box(which = "plot", bty = "l")
mtext('this study', side = 1, line = 2.5, cex = 0.75)
mtext('Irisarri et al. 2017', side = 2, line = 2.5, cex = 0.75)

abline(a = 0, b = 1, lty = 2)

points(branching.times(fulltree_fam)[as.character(xx[,1])], branching.times(irisarri_fam)[as.character(xx[,2])], pch = 21, bg = gray(0.85), cex = 1.5)

correlation <- cor.test(branching.times(fulltree_fam)[as.character(xx[,1])], branching.times(irisarri_fam)[as.character(xx[,2])])

# mtext(bquote(italic("r") == .(round(correlation$estimate, 2))), side = 3, line = 1.2)
text(grconvertX(0.90, from = 'npc', to = 'user'), grconvertY(0.05, from = 'npc', to = 'user'), bquote(italic("r") == .(round(correlation$estimate, 2))), xpd = NA, cex = 1.5)

mtext('Irisarri et al. 2017', side = 3, font = 2)

dev.off()


