# figure showing how the full tree was split into a few large clades for TreePL. Also show nodes that calibrated as secondary calibrations, and also smoothing parameters.

library(ape)

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

alldat <- read.csv('./data/alldat.csv')

fulltreefile <- "./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre"
backboneFile <- './data/2.time_calibration_imputation/MCMCTreeAnalyses/mcmctree_useData2/iterA/FigTree.tre'
mcmcPostFile <- './data/2.time_calibration_imputation/MCMCTreeAnalyses/mcmctree_useData2/iterA/mcmc.txt'


tr <- read.tree(fulltreefile)
tr <- ladderize(tr, right = TRUE)

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
tr <- rotate(tr, nn)

tr <- read.tree(text = write.tree(tr)) # otherwise, node rotation doesn't change structure



backboneTree <- MCMCtreeR::readMCMCtree(backboneFile)
mcmcPost <- data.table::fread(mcmcPostFile, data.table = FALSE)

# match up nodes in mcmcPost to nodes in full tree
mcmcCols <- grep('t_n', colnames(mcmcPost), value = TRUE)
mcmcCorrespondingNodes <- c()
for (i in 1:length(mcmcCols)) {
	if (all(geiger::tips(backboneTree[[1]], as.numeric(gsub('t_n', '', mcmcCols[i]))) %in% tr$tip.label)) {
		mcmcCorrespondingNodes[i] <- getMRCA(tr, geiger::tips(backboneTree[[1]], as.numeric(gsub('t_n', '', mcmcCols[i]))))
	} else {
		mcmcCorrespondingNodes[i] <- NA
	}
}




subtreeTips <- list(

	Gekkota = c('Coleonyx_variegatus', 'Crenadactylus_naso', 'Nephrurus_levis', 'Lialis_burtonis', 'Hemidactylus_frenatus', 'Gonatodes_humeralis', 'Thecadactylus_solimoensis'),
	
	Scincoidea = c('Plestiodon_skiltonianus', 'Zonosaurus_laticaudatus', 'Smaug_mossambicus', 'Xantusia_vigilis'),
	
	Laterata = c('Bipes_biporus', 'Monopeltis_guentheri', 'Dryadosaura_nordestina', 'Ameiva_ameiva', 'Lacerta_viridis', 'Rhineura_floridana', 'Ptychoglossus_brevifrontalis', 'Blanus_strauchi', 'Cadea_blanoides', 'Trogonophis_wiegmanni'),
	
	Iguania = c('Chamaeleo_calyptratus', 'Moloch_horridus', 'Uma_notata', 'Dipsosaurus_dorsalis', 'Anolis_punctatus', 'Diplolaemus_darwinii', 'Polychrus_liogaster', 'Hoplocercus_spinosus', 'Leiocephalus_macropus', 'Phymaturus_palluma', 'Oplurus_cyclurus', 'Basiliscus_basiliscus', 'Crotaphytus_collaris', 'Stenocercus_squarrosus'),
	
	Anguimorpha = c('Heloderma_suspectum', 'Ophiodes_striatus', 'Anniella_pulchra', 'Elgaria_kingii', 'Xenosaurus_platyceps', 'Shinisaurus_crocodilurus', 'Varanus_griseus', 'Lanthanotus_borneensis'),
	
	Serpentes = c('Acrochordus_granulatus', 'Anilius_scytale', 'Liotyphlops_ternetzii', 'Charina_bottae', 'Casarea_dussumieri', 'Contia_tenuis', 'Cylindrophis_ruffus', 'Micrurus_fulvius', 'Gerrhopilus_mirus', 'Fordonia_leucobalia', 'Atractaspis_fallax', 'Namibiana_labialis', 'Loxocemus_bicolor', 'Python_regius', 'Trachyboa_boulengeri', 'Anilios_waitii', 'Rhinophis_saffragamus', 'Agkistrodon_contortrix', 'Xenopeltis_unicolor', 'Anomochilus_leonardi', 'Pareas_nuchalis', 'Kladirostratus_acutus', 'Fimbrios_klossi', 'Xenophidion_schaeferi', 'Xenotyphlops_grandidieri')
	
)

all(unlist(subtreeTips) %in% tr$tip.label)


topSmoothVal <- numeric(length(subtreeTips)+1)
names(topSmoothVal) <- c(names(subtreeTips), 'scaffold')

# hard code best smoother, based on manual examination
topSmoothVal['Gekkota'] <- 1
topSmoothVal['Scincoidea'] <- 1
topSmoothVal['Laterata'] <- 1
topSmoothVal['Iguania'] <- 0.1
topSmoothVal['Anguimorpha'] <- 0.1
topSmoothVal['Serpentes'] <- 0.1
topSmoothVal['scaffold'] <- 1


backboneTree <- backboneTree[[1]]


# identify branches that are consistent with constraint tree
commonsp <- intersect(backboneTree$tip.label, tr$tip.label)

edgeVec <- sapply(tr$edge[,2], function(x) any(commonsp %in% geiger::tips(tr, x)))
edgeVec <- as.character(edgeVec)
table(edgeVec)
edgeVec[which(edgeVec == 'TRUE')] <- 'darkorchid3'
edgeVec[which(edgeVec == 'FALSE')] <- 'gray80'
table(edgeVec)

tipCol <- rep(gray(0.85), Ntip(tr))
tipCol[tr$tip.label %in% backboneTree$tip.label] <- 'black'
table(tipCol)

scaffoldEdges <- which(!tr$edge[,2] %in% unlist(lapply(subtreeTips, function(x) phangorn::Descendants(tr, getMRCA(tr, x), type = 'all'))))

scaffoldEdgeCols <- rep('gray80', Nedge(tr))
scaffoldEdgeCols[scaffoldEdges] <- 'darkorchid3'

edgeWidths <- rep(0.4, Nedge(tr))
edgeWidths[scaffoldEdges] <- 2


# function that returns coordinates that outline a clade in a phylogeny
getCladePoly <- function(tree, node) {
	
	spanningInd <- range(which(tree$tip.label %in% geiger::tips(tree, node)))
	# spanningTips <- getSpanningTips(tree, node)
	
	# spanningInd <- which(tree$tip.label %in% spanningTips)
	
	nn <- min(spanningInd)

	nodePath <- nn
	while (nn != node) {
		nn <- tree$edge[tree$edge[,2] == nn, 1]
		nodePath <- c(nodePath, nn)
	}

	xx1 <- c()
	yy1 <- c()
	for (i in 1:length(nodePath)) {
				
		if (i > 1) {
			xx1 <- c(xx1, get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[nodePath[i]])
			desc <- phangorn::Descendants(tree, nodePath[i], type = 'children')
			yy1 <- c(yy1, min(get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[desc]))
		}
				
		xx1 <- c(xx1, get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[nodePath[i]])
		yy1 <- c(yy1, get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[nodePath[i]])
		
	}
	
	nn <- max(spanningInd)

	nodePath <- nn
	while (nn != node) {
		nn <- tree$edge[tree$edge[,2] == nn, 1]
		nodePath <- c(nodePath, nn)
	}

	xx2 <- c()
	yy2 <- c()
	for (i in 1:length(nodePath)) {
		
		if (i > 1) {
			xx2 <- c(xx2, get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[nodePath[i]])
			desc <- phangorn::Descendants(tree, nodePath[i], type = 'children')
			yy2 <- c(yy2, max(get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[desc]))
		}
		
		xx2 <- c(xx2, get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[nodePath[i]])
		yy2 <- c(yy2, get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[nodePath[i]])
		
	}
	
	xx <- c(xx1, rev(xx2), xx1[1])
	yy <- c(yy1, rev(yy2), yy1[1])
	
	return(cbind(xx, yy))
	
}	
	




png('~/Downloads/figS6.png', width = 7, height = 9, units = 'in', res = 400)

par(mar = c(1,0,1,4))

plot.phylo(tr, show.tip.label = FALSE, edge.col = scaffoldEdgeCols, edge.width = edgeWidths, tip.color = tipCol, cex = 0.5)

segments(x0 = 25, x1 = 50, y0 = 1300, y1 = 1300, lwd = 2, col = 'darkorchid3')
text(x = 37.5, y = 1300, pos = 1, 'scaffold (1)', cex = 0.75)

points(x = 25, y = 800, pch = 21, cex = 0.5, bg = 'dodgerblue', lwd = 0.5)
text(x = 25, y = 800, 'secondary calibration from MCMCtree analysis', pos = 4, cex = 0.75)

for (i in 1:length(subtreeTips)) {
	
	cladeMRCA <- getMRCA(tr, subtreeTips[[i]])
	
	coords <- getCladePoly(tr, cladeMRCA)
	
	polygon(coords, col = adjustcolor('red', alpha.f = 0.05), border = 'red', lwd = 0.8)	
	
	spanningInd <- which(tr$tip.label %in% geiger::tips(tr, cladeMRCA))
	
	lab <- paste0(names(subtreeTips)[i], '\n(', topSmoothVal[names(subtreeTips)[i]], ')')
	text(x = max(branching.times(tr)), y = mean(spanningInd), lab, pos = 4, xpd = NA, cex = 0.75)
	
	
}

nodelabels(node = mcmcCorrespondingNodes, pch = 21, cex = 0.5, bg = 'dodgerblue', lwd = 0.5)


dev.off()




	



