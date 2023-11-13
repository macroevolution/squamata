rm(list = ls())

library(ape)
library(phytools)
library(phangorn)
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

get_discord <- function(pgen, tgen) {
  discord = rep(FALSE, Nnode(pgen))
  for (i in 1:length(discord)) {
    nodenum = i + Ntip(pgen)
    desc = pgen$tip.label[ Descendants(pgen, nodenum, type = "tips")[[1]] ]
    nodenum2 = findMRCA(tgen, desc)
    desc2 = tgen$tip.label[ Descendants(tgen, nodenum2, type = "tips")[[1]] ]
    diff = c(setdiff(desc, desc2), setdiff(desc2, desc))
    if (length(diff) > 0) {
      discord[i] = TRUE
    }
  }
  cat(table(discord))
  return(discord)
}

# our tree
p = read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")

# Tonini
t = read.tree("phylogenetic_inference/treesForComparison/tonini2016.tre")
t = drop.tip(t, "Sphenodon_punctatus")
t = read.tree(text = write.tree(ladderize(t)))
keep = intersect(t$tip.label, p$tip.label)
t1 = keep.tip(t, keep)
p1 = keep.tip(p, keep)

fam = get_subsample(alldat, t1, "family")
tfam = keep.tip(t1, fam)
pfam = keep.tip(p1, fam)
tfam$tip.label = names(fam[ match(tfam$tip.label, fam) ])
pfam$tip.label = names(fam[ match(pfam$tip.label, fam) ])
comp = cophylo(pfam, tfam)
tfam2 = minRotate(tfam, rev(match(pfam$tip.label, tfam$tip.label)))

gen = get_subsample(alldat, t1, "genus")
tgen = keep.tip(t1, gen)
pgen = keep.tip(p1, gen)
dt = get_discord(pgen, tgen)

# Burbrink
b = read.tree("phylogenetic_inference/treesForComparison/Burbrink2020.tre")
b = drop.tip(b, "Sphenodon_punctatus")

# will compare full tree
keep = intersect(b$tip.label, p$tip.label)
b1 = keep.tip(b, keep)
p1 = keep.tip(p, keep)
db = get_discord(p1, b1)

# family level tree
gen = gsub("_\\S+", "", b$tip.label)
fams = alldat[match(gen, alldat$genus), "family"]
names(fams) = b$tip.label
fams = fams[complete.cases(fams)]
fams2 = unique(fams)
keep = rep(NA, length(fams2))
for (i in 1:length(fams2)) {
  keep[i] = sample(names(fams[fams == fams2[i]]), 1)
}
bfam = keep.tip(b, keep)
bfam$tip.label = fams[ match(bfam$tip.label, names(fams)) ]
both = intersect(bfam$tip.label, pfam$tip.label)
bfam2 = keep.tip(bfam, both)
pfam2 = keep.tip(pfam, both)
comp2 = cophylo(pfam2, bfam2)


pdf("manuscript/figures/SupplInfo_FigS39.pdf", width = 10, height = 13)
layout(mat = matrix(c(1, 5, 2, 6, 3, 7, 4, 8), nrow = 2), 
       widths = c(0.4, 0.25, 0.1, 0.25))
###### tonini
par(mar = c(3, 0, 0, 0), xpd = TRUE)
plot(pgen, show.tip.label = F, lwd = 0.01, edge.color = "black")
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
text(labels = "A", x = pp$x.lim[2] * 0.01, y = pp$y.lim[2], font = 2, cex = 2)
nodelabels("", frame = "none", pch = 16, col = "red",
           cex = ifelse(dt == TRUE, 1, 0.01))
axisPhylo()
par(mar = c(0, 0, 0, 7), xpd = TRUE)
xx = comp$assoc
t1 = comp$trees[[1]]
t2 = comp$trees[[2]]
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
par(mar = c(0, 7, 0, 0), xpd = TRUE)
plot.phylo(t2, direction = "leftwards", show.tip.label = F)
tiplabels(t2$tip.label, adj = c(1, 0.5), frame = "none")

######### burbrink
par(mar = c(3, 0, 0, 0), xpd = TRUE)
plot(p1, show.tip.label = F, lwd = 0.01, edge.color = "black")
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
text(labels = "B", x = pp$x.lim[2] * 0.01, y = pp$y.lim[2], font = 2, cex = 2)
nodelabels("", frame = "none", pch = 16, col = "red",
           cex = ifelse(db == TRUE, 1, 0.01))
axisPhylo()
par(mar = c(0, 0, 0, 7), xpd = TRUE)
xx = comp2$assoc
t1 = comp2$trees[[1]]
t2 = comp2$trees[[2]]
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
par(mar = c(0, 7, 0, 0), xpd = TRUE)
plot.phylo(t2, direction = "leftwards", show.tip.label = F)
tiplabels(t2$tip.label, adj = c(1, 0.5), frame = "none")

dev.off()

