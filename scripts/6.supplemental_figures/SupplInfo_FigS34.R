library(ape)
library(dplyr)
library(tidyr)

setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

t = read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
d = read.csv("Data_S2.csv")
mt = c("COI", "ND1", "ND4", "rRNA_12S", "rRNA_16S", "CYTB", "ND2")

loc = unique(d$locus)
x = matrix(NA, nrow = length(t$tip.label), 
           ncol = length(loc), dimnames = list(t$tip.label, loc))
for (i in 1:nrow(x)) {
 d1 = d[which(d$finalTipName == rownames(x)[i]), ]
 d2 = d1[match(loc, d1$locus), "accession"]
 row = d2
 row[complete.cases(d2)] = 1
 row[grepl("SqCL", d2)] = 3
 row[loc %in% mt & complete.cases(d2)] = 2
 row[is.na(d2)] = 0
 x[i, ] = row
}

m = apply(x, 2, as.numeric)
counts = apply(m, 2, function(x) (sum(x > 0)))
m2 = m[,names(counts)[order(counts)]]
m2 = t(m2)

alldat <- read.csv('./data/alldat.csv')
edgecols = rep("gray70", nrow(t$edge))
gen = alldat[alldat$inTopoConstraint == TRUE, "treename"]
for (i in 1:length(gen)) {
  tip = which(t$tip.label == gen[i])
  edgecols[ which(t$edge[, 2] == tip) ] = "#fc8d62"
}

pdf("~/Desktop/supermatrix.pdf", width=6, height=8)
layout(matrix(c(1, 2), ncol = 2), widths = c(0.35, 0.65))
par(mar = c(2.6, 0, 0, 1))
plot.phylo(t, direction = "rightwards", show.tip.label = FALSE, 
           edge.width = 0.2, edge.color = edgecols)
lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
par(mar = c(4, 0, 1.3, 1))
gl <- 1:(nrow(m2) + 1)
image(1:dim(m2)[1], 1:dim(m2)[2], m2, 
      col = c("gray90", "#66c2a5", "#8da0cb", "#fc8d62"), axes=F, xlab="")
xvals = pretty(c(1,dim(m2)[1]))[1:6]
axis(1, at=xvals, labels=xvals)
mtext("number of loci", 1, line=2.5)
dev.off()
