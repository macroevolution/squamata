rm(list = ls())

library(ape)
library(phangorn)
library(phytools)
setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

alldat <- read.csv('./data/alldat.csv')

get_subsample <- function(d, t, type) {
  # first get the list of genera
  gen = unique(d[ , type])
  # randomly sample one tip per fam
  keep = rep(NA, length(gen))
  
  for (i in 1:length(gen)) {
    sps = d[d[, type] == gen[i], "treename"]
    sps = sps[sps %in% t$tip.label]
    if (length(sps) > 0) {
      keep[i] = sample(sps, 1)
    }
  }
  names(keep) = gen
  keep2 = keep[complete.cases(keep)]
  
  return(keep2)
}



con = read.tree("./data/1.tree_inference/mainTrees/fulltree_default_con_1_raxmlOpt.raxml.final.tre")
uncon = read.tree("./data/1.tree_inference/mainTrees/fulltree_default_uncon_2_raxmlOpt.raxml.bestTree")

con = root(con, "Sphenodon_punctatus", resolve.root = T)
con = drop.tip(con, "Sphenodon_punctatus")
uncon = root(uncon, "Sphenodon_punctatus", resolve.root = T)
uncon = drop.tip(uncon, "Sphenodon_punctatus")
con = read.tree(text = write.tree(ladderize(con)))
uncon = read.tree(text = write.tree(ladderize(uncon)))

# constrained, but ultrametric
con2 = read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")

# family or genus tree with node circles
gen = get_subsample(alldat, con2, "genus")
gencon = keep.tip(con2, gen)
genuncon = keep.tip(uncon, gen)
discord = rep(FALSE, Nnode(gencon))
for (i in 1:length(discord)) {
  nodenum = i + Ntip(gencon)
  desc = gencon$tip.label[ Descendants(gencon, nodenum, type = "tips")[[1]] ]
  nodenum2 = findMRCA(genuncon, desc)
  desc2 = genuncon$tip.label[ Descendants(genuncon, nodenum2, type = "tips")[[1]] ]
  diff = c(setdiff(desc, desc2), setdiff(desc2, desc))
  if (length(diff) > 0) {
    discord[i] = TRUE
  }
}
table(discord) / Nnode(gencon)


# family cophylo
fam = get_subsample(alldat, con, "family")
famcon = keep.tip(con, fam)
famuncon = keep.tip(uncon, fam)
famcon$tip.label = names(fam[ match(famcon$tip.label, fam) ])
famuncon$tip.label = names(fam[ match(famuncon$tip.label, fam) ])
comp = cophylo(famcon, famuncon)

# can get rotated tips as comp$trees[[1]] & comp$trees[[2]]
t1 = comp$trees[[1]]
t2 = comp$trees[[2]]
xx = comp$assoc

pdf("manuscript/figures/SupplInfo_FigS35.pdf", width = 10, height = 8)
layout(mat = matrix(c(1, 2, 3, 4), nrow = 1), 
       widths = c(0.5, 0.2, 0.1, 0.2))
par(mar = c(3, 0, 0, 0), xpd = TRUE)
plot(gencon, show.tip.label = F, lwd = 0.01, edge.color = "black")
nodelabels("", frame = "none", pch = 16, col = "red",
           cex = ifelse(discord == TRUE, 1, 0.01))
axisPhylo()
par(mar = c(0, 0, 0, 6), xpd = TRUE)
plot(t1, show.tip.label = F)
tiplabels(t1$tip.label, adj = c(0, 0.5), frame = "none")
par(mar = c(0, 0, 0, 0), xpd = TRUE)
plot(NA, xlim = c(0, 1), ylim = c(1, Ntip(t1)), axes = F)
for (i in 1:nrow(xx)) {
 y1 = which(t1$tip.label == xx[i, 1]) 
 y2 = which(t2$tip.label == xx[i, 1])
 if (!identical(y1, integer(0)) && !identical(y2, integer(0)) ) {
   if (y1 == y2) {
     color = "black"
     alpha = 0.5
   } else {
     color = "red"
     alpha = 1
   }
   lines(x = c(0.1, 0.9), y = c(y1, y2), col = alpha(color, alpha))
 }
}
par(mar = c(0, 6, 0, 0), xpd = TRUE)
plot.phylo(t2, direction = "leftwards", show.tip.label = F)
tiplabels(t2$tip.label, adj = c(1, 0.5), frame = "none")
dev.off()

