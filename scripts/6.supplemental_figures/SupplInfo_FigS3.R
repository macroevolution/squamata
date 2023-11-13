rm(list = ls())
setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

library(ape)
library(scales)

t = read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
d = read.csv("./data/squamatesCoreTaxonomy_treeTaxa.csv")

edgecols = rep(alpha("gray50", 0.01), nrow(t$edge))
gen = d[d$inTopoConstraint == TRUE, "treename"]
for (i in 1:length(gen)) {
  tip = which(t$tip.label == gen[i])
  anc = c(tip, phangorn::Ancestors(t, tip, type = "all"))
  edgecols[ which(t$edge[, 2] %in% anc) ] = "coral1"
}

pdf("manuscript/figures/SupplInfo_FigS3.pdf", height = 10, width = 10)
par(mar = c(0, 0, 0, 0), xpd = T)
plot(t, show.tip.label = F, edge.color = "gray70",
     lwd = 0.1, type = "fan")
par(new=TRUE)
plot.phylo(t, show.tip.label = F, edge.color = edgecols,
     lwd = 1, type = "fan")
dev.off()