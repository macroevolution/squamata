require(BAMMtools)

setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

alldat <- read.csv('./data/alldat.csv')

tr <- read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
tr <- ladderize(tr)

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(tr, intersect(tr$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
tr <- rotate(tr, nn)

tr <- read.tree(text = write.tree(tr)) # otherwise, node rotation doesn't change structure



snakes <- alldat[alldat$clade_Serpentes == 1, 'treename']
lizards <- alldat[alldat$clade_Serpentes == 0, 'treename']


# Proportionally thinned version of tree, cohort analysis with notable clades highlighted

ed2 <- getEventData(tr, eventdata = './data/3.trait_preparation_diversification/bestBAMMrun/event_data.txt', nsamples = 1000, burnin = 0.1)


nPreserved <- 150

preservedTips <- seq(1, Ntip(tr), by = round(Ntip(tr) / nPreserved))
preservedTips <- rev(tr$tip.label)[preservedTips] # rev to retain Dibamidae
droppedTips <- setdiff(tr$tip.label, preservedTips)

subtree <- keep.tip(tr, preservedTips)

subtree <- ladderize(subtree, right = TRUE)

# rotate Iguania node to avoid Liolaemus being next to snakes and giving the allusion that rates are high for a clade
nn <- getMRCA(subtree, intersect(subtree$tip.label, alldat[alldat[, 'clade_Iguania'] == 1, 'treename']))
subtree <- rotate(subtree, nn)

subtree <- read.tree(text = write.tree(subtree)) # otherwise, node rotation doesn't change structure

tipList <- lapply(preservedTips, function(x) x)
names(tipList) <- preservedTips
patdist <- cophenetic.phylo(tr)
tol <- 1e-6
for (i in 1:length(droppedTips)) {
	xx <- patdist[droppedTips[i], preservedTips]
	assoc <- names(xx)[which(abs(xx - min(xx)) < tol)]
	# if (length(assoc) > 1) stop()
	for (j in 1:length(assoc)) {
		tipList[[assoc[j]]] <- c(tipList[[assoc[j]]], droppedTips[i])
	}
}


subtreeBamm <- subtreeBAMM(ed2, tips = subtree$tip.label)
subtreeBamm <- dtRates(subtreeBamm, tau = 0.001)

# cohort matrix
xx <- getCohortMatrix(subtreeBamm)

png('~/Downloads/figS26.png', width = 7, height = 7, units = 'in', res = 500)

x <- xx
ephy <- subtreeBamm
ofs = -0.02
figs <- matrix(c(0,0.2,0.8,1,
            	 0.2,0.95,0.8+ofs,1,
                 0,0.2-ofs,0,0.8,
                 0.2,0.95,0,0.8,
                 0.95,1,0.25,0.75
                 ), byrow=TRUE,
               nrow=5, ncol=4);
useraster <- ifelse(dim(x)[1] > 500, TRUE, FALSE);
pal <- "RdYlBu";
col <- colorRampPalette(BAMMtools:::palettes[["RdYlBu"]])(64)
# col <- colorRampPalette(c('black', 'white'))(64)
col <- colorRampPalette(c('white', 'black'))(64)
# col <- viridis::rocket(64)
ncolors <- length(col);
breaks <- quantile(seq(0,1.01,length.out=100),probs=seq(0,1,length.out=ncolors+1));
	
index <- match(ephy$tip.label, rownames(x));
x <- x[index, index];



# label major nodes
clades <- c('Serpentes', 'Alethinophidia', 'Colubriformes', 'Iguania', 'Gekkonidae', 'Scincidae', 'Liolaemus')

nn <- numeric(length(clades))
for (i in 1:length(clades)) {
	
	if (paste0('clade_', clades[i]) %in% colnames(alldat)) {
		sp <- alldat[alldat[, paste0('clade_', clades[i])] == 1, 'treename']
	} else if (clades[i] == 'lizards') {
		sp <- alldat[alldat[, 'clade_Serpentes'] == 0, 'treename']
	} else if (clades[i] %in% alldat$family) {
		sp <- alldat[alldat$family == clades[i], 'treename']
	} else if (clades[i] %in% alldat$genus) {
		sp <- alldat[alldat$genus == clades[i], 'treename']
	} else if (clades[i] == 'Oz_Sphenomorphines') {
		sp <- alldat[which(as.character(alldat$geogRadiation) == 'spheno-oz'), 'treename']
	} else if (clades[i] == 'limb_reduced_lizards') {
		sp <- alldat[alldat$numberLimbs < 4 | alldat$numberDigits < 10, 'treename']
		sp <- setdiff(sp, alldat[alldat[, 'clade_Serpentes'] == 1, 'treename'])
	}
		
	nn[i] <- getMRCA(subtree, intersect(sp, subtree$tip.label))
}


# left side tree
par(fig = figs[3,], new=FALSE, mar = c(5,1,0,0));
plot(ephy, pal = pal, lwd = 1.3, direction="rightwards", breaksmethod = 'jenks')

nodelabels(1:length(clades), node = nn, frame = 'circle', cex = 0.6, bg = adjustcolor('light blue', alpha.f = 0.65), font = 2)

# top tree
par(fig = figs[2,], new = TRUE, mar = c(0,0,1,4));
phy <- BAMMtools:::as.phylo.bammdata(ephy);
bt <- max(ephy$end)
plot.phylo(phy, edge.width = 1.3, direction="downwards", show.tip.label = FALSE, x.lim = c(1,length(phy$tip.label)), y.lim = c(0,bt), edge.color = gray(0.35))


nodelabels(1:length(clades), node = nn, frame = 'circle', cex = 0.6, bg = adjustcolor('light blue', alpha.f = 0.65), font = 2)

text(x = grconvertX(0.85, from = 'nfc', to = 'user'), y = grconvertY(seq(0.85, 0.35, length.out = length(clades)), from = 'nfc', to = 'user'), labels = paste0(1:length(clades), '. ', clades), cex = 0.75, pos = 4, xpd = NA)

par(fig = figs[4,], new=TRUE, mar = c(5,0,0,4));
plot(0,0,type="n",axes=FALSE,ann=FALSE,xlim=c(0,1),ylim=c(0,1))

image(x,axes=FALSE,xlab="",ylab="",col=col,xlim=c(0,1),ylim=c(0,1),breaks=breaks,add=TRUE,useRaster=useraster);

BAMMtools:::barLegend(col, breaks,fig=c(0.90, 0.93, 0.35, 0.65), side = 2, cex.axis = 0.75);

text(x = grconvertX(0.93, from = 'ndc', to = 'user'), y = grconvertY(0.5, from = 'ndc', to = 'user'), 'pairwise prob.\nof shared rate regime', srt = -90, cex = 0.7, xpd = NA, pos = 3)

dev.off()

