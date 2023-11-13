rm(list = ls())
library(ape)

read_tree <- function(x) {
  t = read.tree(x)

  # root the tree
  outs = c("hg38", "galGal5", "taeGut2", "chrPic1", "allMis1",
           "32244", "ISIS373002", "UMFS-10956c", "H20145a",
           "HCD-2620a", "alligator_mississippiensis",
           "H2662a", "chrysemys_picta", "LM-67b",
           "gallus_gallus", "taeGut1", "galGal3", "hg19")
  outs2 = outs[outs %in% t$tip.label]
  t2 = root(t, outs2[1], resolve.root = T)
  t3 = drop.tip(t2, outs2)
  
  # map the tree to species names
  taxonTable <- read.csv('./data/0.phylogenomic_data/sqCL_taxa_repDBmatching.csv', stringsAsFactors=FALSE)
  ## manual matching supercedes repDB auto matching, unless no manual match found
  taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'repDBtaxon']
  head(taxonTable)
  table(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch))
  ## for the few individuals that don't have a good match
  ## just match it to the existing species name
  taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'manualMatch'] <- taxonTable[which(taxonTable$manualMatch == '' | is.na(taxonTable$manualMatch)), 'species']
  
  miss = t3$tip.label[! (t3$tip.label %in% taxonTable$sample) ]
  t4 = drop.tip(t3, miss)
  t4$tip.label = gsub(" ", "_", taxonTable[ match(t4$tip.label, taxonTable$sample), "manualMatch"])

  t5 = read.tree(text = write.tree(ladderize(t4)))
  return(t5)
}

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

get_discord <- function(t1, t2) {
  discord = rep(FALSE, Nnode(t1))
  for (i in 1:length(discord)) {
    nodenum = i + Ntip(t1)
    desc = t1$tip.label[ phangorn::Descendants(t1, nodenum, type = "tips")[[1]] ]
    nodenum2 = phytools::findMRCA(t2, desc)
    desc2 = t2$tip.label[ phangorn::Descendants(t2, nodenum2, type = "tips")[[1]] ]
    diff = c(setdiff(desc, desc2), setdiff(desc2, desc))
    if (length(diff) > 0) {
      discord[i] = TRUE
    }
  }
  return(discord)
}

setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

# get all trees & root & rename tips
et = read_tree("./data/0.phylogenomic_data/all.ind0.01_loci0.05_n5185.examl.tre")
it = read_tree("./data/0.phylogenomic_data/all.concat_ind0.01_loci0.05_n5185.iqtree.tre")
at = read_tree("./data/0.phylogenomic_data/all_0.01_0.05_10_100_5e-05_SH.astral.tre")

# subsample to genus
alldat <- read.csv('./data/alldat.csv')
gen = get_subsample(alldat, it, "subfamily")
et2 = keep.tip(et, gen)
it2 = keep.tip(it, gen)
at2 = keep.tip(at, gen)

et2$tip.label = names(gen[ match(et2$tip.label, gen) ])
it2$tip.label = names(gen[ match(it2$tip.label, gen) ])
at2$tip.label = names(gen[ match(at2$tip.label, gen) ])
comp1 = phytools::cophylo(it2, at2)
comp2 = phytools::cophylo(it2, et2)

pdf("manuscript/figures/SupplInfo_FigS31.pdf", width = 8, height = 16)
layout(mat = matrix(c(1, 2, 3,
                      4, 5, 6), nrow = 2, byrow = T), 
       widths = c(0.4, 0.2, 0.4))
par(mar = c(0, 0, 0, 6), xpd = TRUE)

xx = comp1$assoc
# can get rotated tips as comp$trees[[1]] & comp$trees[[2]]
t1 = comp1$trees[[1]]
t2 = comp1$trees[[2]]
d1 = get_discord(t2, t1)


plot(t1, show.tip.label = F)
mtext("A", side = 3, line = -1, at = c(-1, 0))
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
    lines(x = c(0.1, 0.9), y = c(y1, y2), col = adjustcolor(color, alpha.f = alpha))
  }
}
par(mar = c(0, 6, 0, 0), xpd = TRUE)
plot.phylo(t2, direction = "leftwards", show.tip.label = F)
tiplabels(t2$tip.label, adj = c(1, 0.5), frame = "none")
nodelabels("", frame = "none", pch = 21, bg = "red",
           cex = ifelse(d1 == TRUE, 1, 0.01))


xx = comp2$assoc
# can get rotated tips as comp$trees[[1]] & comp$trees[[2]]
t1 = comp2$trees[[1]]
t2 = comp2$trees[[2]]
d1 = get_discord(t2, t1)

par(mar = c(0, 0, 0, 6), xpd = TRUE)
plot(t1, show.tip.label = F)
mtext("B", side = 3, line = -1, at = c(-1, 0))
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
    lines(x = c(0.1, 0.9), y = c(y1, y2), col = adjustcolor(color, alpha.f = alpha))
  }
}
par(mar = c(0, 6, 0, 0), xpd = TRUE)
plot.phylo(t2, direction = "leftwards", show.tip.label = F)
tiplabels(t2$tip.label, adj = c(1, 0.5), frame = "none")
nodelabels("", frame = "none", pch = 21, bg = "red",
           cex = ifelse(d1 == TRUE, 1, 0.01))

dev.off()
