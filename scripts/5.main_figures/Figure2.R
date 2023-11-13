#############################################
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, USA.
### July 2023.
#################
### Script goals:
### To plot major phase shifts in trait evolutionary dynamics along squamate evolution.

## 1. Setting things up ----

## Packages:
library(ape)

## Dir:
setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

## Tree to use:
v <- read.tree("./data/5.main_figures/fig2/subfamily-pairs-squams-r1.tre")

## Data:
xx <- read.csv("./data/5.main_figures/fig2/modelFit_cpo_data.csv")
xx$trait <- R.utils::capitalize(xx$trait) ## Capitalize traits (for plots).

## Clade info:
snakes <- extract.clade(v, node = getMRCA(v, tip = c("Grayia_smithii", "Anilios_howi")))$tip.label
lizards <- setdiff(v$tip.label, snakes)
aleth <- extract.clade(v, node = getMRCA(v, tip = c("Anilius_scytale", "Grayia_smithii")))$tip.label
afrop <- extract.clade(v, node = getMRCA(v, tip = c("Rhinophis_phillipsi", "Grayia_smithii")))$tip.label
xx$isSnake <- xx$spanningTaxon1 %in% snakes & xx$spanningTaxon2 %in% snakes
xx$isAlethinophidia <- xx$spanningTaxon1 %in% aleth & xx$spanningTaxon2 %in% aleth
xx$isAfrophidia <- xx$spanningTaxon1 %in% afrop & xx$spanningTaxon2 %in% afrop

## Nodes of interest to highlight with different colors in plots:
xx$serpentes <- 0
xx$alethinophidia <- 0
xx$serpentes[xx$spanningTaxon1 == "Gyalopion_quadrangulare" & xx$spanningTaxon2 == "Anilios_howi"] <- 1
xx$alethinophidia[xx$spanningTaxon1 == "Gyalopion_quadrangulare" & xx$spanningTaxon2 == "Tropidophis_melanurus"] <- 1
xx$alethinophidia[xx$spanningTaxon1 == "Gyalopion_quadrangulare" & xx$spanningTaxon2 == "Chilabothrus_strigilatus"] <- 1

## Tips in tree:
tipset <- c(xx$spanningTaxon1, xx$spanningTaxon2)
tipset <- unique(tipset[!is.na(tipset)])

## 2. AIC and r2 improvements by shift ----

## Data:
dsets <- unique(xx$trait)

## Two alternative approaches:
mod <- xx[xx$analysis == "modelFitting" & !is.na(xx$r2diffs), ] 
cpo <- xx[xx$analysis == "CPO" & !is.na(xx$r2diffs), ]

## New dataframe only with the data we care about:
rr <- mod[1:5, ]
rrcpo <- cpo[1:5, ]
for (ii in 2:length(dsets)){
  tmp <- mod[mod$trait == dsets[ii], ]
	rr <- rbind(rr, tmp[1:5, ]) 	## "model fitting" approach.
	tmp <- cpo[cpo$trait == dsets[ii], ]
	rrcpo <- rbind(rrcpo, tmp[1:5, ]) ## CPO approach.
}

## 3. Set up plots ----

## Prune tree to keep only those tips in top five shifts per clade:
tipset <- c(rr$spanningTaxon1, rr$spanningTaxon2)
tipset <- unique(tipset[!is.na(tipset)])

## Keep subset of tips not represented by shifts:
tipdiff <- setdiff(v$tip.label, tipset)
set.seed(123)
tipdiff <- sample(tipdiff, size = 95, replace = FALSE)

## Tree to plot:
v2 <- keep.tip(v, c(tipdiff, tipset))

## Focal traits to plot:
traits <- c("Diet", "Elongation", "Skull", "Mass", "Climate")

## Function: Setting up plots for each trait:
plot.setup <- function(miny, maxy) {
	plot.new()
	par(mgp = c(3, 0.5, 0), mar = c(2, 1, 1, 1)) ## Bottom, left, top, right.
	plot.window(xlim = c(0, 5.5), ylim = c(miny, maxy))
	axis(1, at = seq(-1, 6, by = 1), labels = TRUE, cex.axis = 0.5, tck = -0.05)
	axis(2, at = seq(-0.2, 1, by = 0.2), las = 1 , cex.axis = 0.5, tck = -0.05)
}

## Function: Colors with transparency (for AIC on nodes):
transparentColor <- function(namedColor, alpha = 0.7) {
	res <- col2rgb(namedColor) / 255;
	return(rgb(red = res[1, ], green = res[2, ], blue = res[3, ], alpha = alpha))
}
snake_col <- transparentColor("red")
aleth_col <- transparentColor("orange")
other_col <- transparentColor("gray20")

## Colors of focal clades in plots:
rr$pcol <- rep(other_col, nrow(rr))
rr$pcol[rr$serpentes == 1] <- snake_col
rr$pcol[rr$alethinophidia == 1] <- aleth_col
rrcpo$pcol <- rep(other_col, nrow(rrcpo))
rrcpo$pcol[rrcpo$serpentes == 1] <- snake_col
rrcpo$pcol[rrcpo$alethinophidia == 1] <- aleth_col

## 4. Plot Figure 2 ----

## File:
pdf("figure_parts/v2/Fig_2_raw.pdf", width = 7, height = 2)

## Horizontal layout, row widths different:
lmat <- rbind(1:5, 1:5, 6:10, 6:10, 6:10)
layout(lmat)
par(mar = c(0.5, 0.5, 0.5, 0.5), 
    oma = c(0, 0.75, 1, 0)) ## Bottom, left, top, right.

## Plot tree, looping over traits:
for (ii in 1:length(traits)) {
 	tset <- rr[rr$trait == traits[ii], ] ## Data for focal trait only.
 	tset$normAIC <- tset$aicdiffs / max(tset$aicdiffs) ## Normalizing AICs by the max of each trait.
	
	## Plot tree:
 	plot.phylo(ladderize(v2, right = FALSE), 
 	           show.tip.label = FALSE, edge.color = "gray70",
 	           edge.width = 0.5,
 	           direction = "upwards")
 	
 	## Plot title (== trait):
	mtext(text = traits[ii], side = 3, line = 0, cex = 0.5, adj = 0.5)
	
	## Plot node annotations, looping over nodes:
	for (kk in 1:nrow(tset)) {
	  nd <- getMRCA(v2, tip = unlist(tset[kk, c("spanningTaxon1", "spanningTaxon2")]))
		if (tset$serpentes[kk] == 1) { bgcol <- snake_col
		} else if (tset$alethinophidia[kk] == 1) { bgcol <- aleth_col 
		} else { bgcol <- other_col	}
		sf <- 2.5 ## Base size of node circles.
		nodelabels(node = nd, pch = 21, bg = bgcol, col = "black", cex = 0.5+(sf*tset$normAIC[kk]), lwd = 0.5)
	}
}

## Axis limits:
minyl <- c(0.0, 0.0, 0.0, 0.0, 0.0)
maxyl <- rep(1, 5)

## Plot r2 scatter plots:
for (ii in 1:length(traits)) {
  
  ## Add points for the "model fitting" approach:
	plot.setup(minyl[ii], maxyl[ii])
  tseta <- rr[rr$trait == traits[ii], ]
	for (kk in 1:nrow(tseta)) {
    points(kk, tseta$r2[kk], pch = 24, cex = 1, lwd = 0.5,
           bg = tseta$pcol[kk], 
           col = "black")           
	}
	## Lines:
	lines(c(0, 1:nrow(tseta)), c(0, tseta$r2), type = "b", 
	      pch = "", lwd = 0.5, lty = 1, 
	      col = "gray30")
	
 	## Add points for the CPO approach:
 	tsetb <- rrcpo[rrcpo$trait == traits[ii], ]
 	for (kk in 1:nrow(tsetb)) {
  		points(kk, tsetb$r2[kk], pch = 22, cex = 1, lwd = 0.5,
  		       bg = tsetb$pcol[kk],
  		       col = "black")
	}
 	## Lines:
 	lines(c(0, 1:nrow(tsetb)), c(0, tsetb$r2), type = "b", pch = "", 
 	      lwd = 0.5, lty = 2, 
 	      col = "gray30")
}
dev.off()

## End of script.
