#setwd('~/Dropbox (Personal)/research_projects-ACTIVE/squamate_tree/')

# working dir is directory with model fitting files and scripts
setwd('~/Dropbox/Oz_Crown_Ages/dataArchive/scripts/4.comparative_analyses/model-fits')

library(ape)
library(RColorBrewer)
library(cutphylo)

source('1.load-data.R')

alldat <- read.csv('../../../data/alldat.csv')

#sfam_tree <- ape::read.tree("scripts/fit-shifts/subfamily-pairs-squams-r1.tre")
#clademat  <- buildSpanningSet_DF(sfam_tree)

#fulltree <- read.tree('./final-trees/best_ultrametric_fulltree_ddBD_revision.tre')
#fulltree <- read.tree(text = write.tree(ladderize(fulltree)))

vals = c("ancDistDiet", "dietRate", 
         "ancDistVert", "vertebralRate",
         "skullAncDist", "skullRate",
         "ancDistElongationIndex", "elongationRate",
         "mass", "ancDistMass", "massRate",
         "climateAncDist", "climateNicheRate")
logit = c(FALSE, TRUE, 
          FALSE, TRUE, 
          FALSE, TRUE,
          FALSE, TRUE,
          TRUE, TRUE, TRUE,
          FALSE, TRUE)

datList2 = vector("list", length(vals))
for (i in 1:length(vals)) {
  val1 = alldat[, vals[i]]
  if (logit[i]) {
    val1 = log(val1)
  }
  names(val1) = alldat$treename
  val1 = val1[complete.cases(val1)]
  
  common <- intersect(names(val1), fulltree$tip.label)
  datList2[[i]] <- as.matrix(val1[common])
}
names(datList2) = vals
datList = c(datList, datList2)

datList <- datList2

order =  c("diet", "ancDistDiet", "dietRate", 
           "vertebral count", "ancDistVert", "vertebralRate",
           "skull shape", "skullAncDist", "skullRate",
           "elongation", "ancDistElongationIndex", "elongationRate",
           "mass", "ancDistMass", "massRate",
           "climate", "climateAncDist", "climateNicheRate")
datList = datList[order]
xlabs = c("multivariate diet", expression(Psi[diet]), expression(TR[diet]), 
          "vertebral count", expression(Psi[vert]), expression(TR[vert]),
          "skull shape", expression(Psi[skull]), expression(TR[skull]),
          "elongation", expression(Psi[elong]), expression(TR[elong]),
          "mass", expression(Psi[mass]), expression(TR[mass]),
          "multivariate climate", expression(Psi[climate]), expression(TR[climate]))
xlabs = c("A. multivariate diet", 
          bquote("B."~Psi[diet]),
          bquote("C."~TR[diet]),
          "D. vertebral count", 
          bquote("E."~Psi[vert]),
          bquote("F."~TR[vert]),
          "G. skull shape", 
          bquote("H."~Psi[skull]),
          bquote("I."~TR[skull]),
          "J. elongation", 
          bquote("K."~Psi[elong]),
          bquote("L."~TR[elong]),
          "M. mass", 
          bquote("N."~Psi[mass]),
          bquote("O."~TR[mass]),
          "P. multivariate climate", 
          bquote("Q."~Psi[climate]),
          bquote("R."~TR[climate]))

nodes =  vector("list", length(datList))
noderes = vector("list", length(datList))
for (i in 1:length(datList)) {
  print(i)
  subtree = keep.tip(fulltree, rownames(datList[[i]]))
  nodes[[i]] = find_clades_pairs(phy = subtree, 
                         fulltree = fulltree, clademat,
                         clip.NA = F) - Ntip(subtree)

  m3 = datList[[i]]
  allmean = colMeans(m3)
  SSE = rep(NA, ncol(m3))
  for (x in 1:ncol(m3)) {
    SSE[x] = sum(sapply(m3[,x], function(y) (y - allmean[x])^2))
  }
  allSSE = sum(SSE)
  
  node_r2 = rep(NA, Nnode(subtree))
  for (j in 1:Nnode(subtree)) {
    tips1 = subtree$tip.label[ phangorn::Descendants(subtree, Ntip(subtree) + j, "tips")[[1]] ]
    tips2 = subtree$tip.label[!(subtree$tip.label %in% tips1)]
    groupmean1 = colMeans(m3[tips1,,drop=FALSE])
    groupmean2 = colMeans(m3[tips2,,drop=FALSE])
    r2_1 = length(tips1) * sum((groupmean1 - allmean)^2) / allSSE
    r2_2 = length(tips2) * sum((groupmean2 - allmean)^2) / allSSE
    node_r2[j] = sum(c(r2_1, r2_2), na.rm = T)
  }
  noderes[[i]] = node_r2
  
  z1 = cut.phylo(subtree, m3, 2, method = "arithmetic")
  print(noderes[[i]][ z1$cuts[2] - Ntip(subtree) ])
  print(z1$r.squared)
  print("****")
}

nodemat = matrix(NA, ncol = length(noderes), nrow = Nnode(sfam_tree))
for (i in 1:ncol(nodemat)) {
  nodemat[ ,i] = c(0, noderes[[i]][ nodes[[i]] ])
}

# remove nodes with NA across all traits
# could also remove nodes with a certain level of missingness (but haven't done it)
sfam_tree$node.label = apply(nodemat, 1, mean, na.rm = T)
drop = unlist(lapply(which(is.na(sfam_tree$node.label)) + Ntip(sfam_tree), function(x) { phangorn::Descendants(sfam_tree, x, type = "tips")[[1]]}))
sfam_tree2 = drop.tip(sfam_tree, sfam_tree$tip.label[drop])
sfam_tree2$node.label[is.na(sfam_tree2$node.label)] = 0

edge.length = rep(0.01, nrow(sfam_tree2$edge))
for (i in 1:Nnode(sfam_tree2)) {
  node = i + Ntip(sfam_tree2)
  edge.length[which(sfam_tree2$edge[,1] %in% node)] = sfam_tree2$node.label[i] 
}
sfam_tree_scale = sfam_tree2
sfam_tree_scale$edge.length = edge.length

snakeNode <- intersect(sfam_tree_scale$tip.label, alldat[alldat$clade_Serpentes == 1, 'treename'])
snakeNode <- getMRCA(sfam_tree_scale, snakeNode)

ColubriformesNode <- intersect(sfam_tree_scale$tip.label, alldat[alldat$clade_Colubriformes == 1, 'treename'])
ColubriformesNode <- getMRCA(sfam_tree_scale, ColubriformesNode)

library(scales)
edge.color = rep(alpha("gray", 0.6), nrow(sfam_tree_scale$edge))
edge.color[which(sfam_tree_scale$edge[,1] %in% c(snakeNode, phangorn::Descendants(sfam_tree_scale, snakeNode, type = 'all')))] = 'darkorange1'
edge.color[which(sfam_tree_scale$edge[,1] %in% c(ColubriformesNode, phangorn::Descendants(sfam_tree_scale, ColubriformesNode, type = 'all')))] = '#386cb0'

ttree = vector("list", length(noderes))
xmax = rep(NA, length(noderes))
for (i in 1:length(noderes)) {
  ttree2 = keep.tip(fulltree, rownames(datList[[i]]))
  ttree2$node.label = noderes[[i]]
  ttree2$node.label[is.na(ttree2$node.label)] = 0
  
  edge.length = rep(0.01, nrow(ttree2$edge))
  for (j in 1:Nnode(ttree2)) {
    node = j + Ntip(ttree2)
    edge.length[which(ttree2$edge[,1] %in% node)] = ttree2$node.label[j] 
  }
  ttree3 = ttree2
  ttree3$edge.length = edge.length
  xmax[i] = max(phytools::nodeHeights(ttree3))
  ttree[[i]] = ttree3
}
xmax = max(xmax)

pdf('manuscript/figures/SupplInfo_FigS15.pdf', width = 6, height = 7)
layout(matrix(seq(1, 18), byrow = TRUE, ncol = 3))
for (i in 1:length(noderes)) {
  ttree3 = ttree[[i]]
  snakeNode <- intersect(ttree3$tip.label, alldat[alldat$clade_Serpentes == 1, 'treename'])
  snakeNode <- getMRCA(ttree3, snakeNode)
  
  ColubriformesNode <- intersect(ttree3$tip.label, alldat[alldat$clade_Colubriformes == 1, 'treename'])
  ColubriformesNode <- getMRCA(ttree3, ColubriformesNode)
  
  tedge.color = rep(alpha("gray", 0.6), nrow(ttree3$edge))
  tedge.color[which(ttree3$edge[,1] %in% c(snakeNode, phangorn::Descendants(ttree3, snakeNode, type = 'all')))] = 'darkorange1'
  tedge.color[which(ttree3$edge[,1] %in% c(ColubriformesNode, phangorn::Descendants(ttree3, ColubriformesNode, type = 'all')))] = '#386cb0'
  par(mar = c(0.5, 0, 2, 0))
  plot.phylo(ttree3, show.node.label = F, show.tip.label = F,
             edge.color = tedge.color, x.lim = c(0, xmax))
  # mtext(text = expression(), side = 3, line = 0, cex = 0.8)
  mtext(xlabs[[i]], side = 3, line = 0, cex = 0.8)
}
dev.off()
