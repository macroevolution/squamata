rm(list = ls())
setwd("~/Dropbox (Personal)/research_projects-ACTIVE/squamate_tree/")

library(ape)

t = read.tree("phylogenetic_inference/data/all_bootstraps/all_IQtree.treefile")
outs = c("chrPic1", "H20145a", "allMis1", "taeGut2",
         "galGal5", "UMFS-10956c", "ISIS373002",
         "hg38")
t1 = drop.tip(t, outs)
t1 = read.tree(text = write.tree(ladderize(t1)))

t1$node.label = as.numeric(t1$node.label)
png("manuscript/figures/SupplInfo_FigS32.png", height = 10, width = 8,
    units = "in", res = 200)
par(mar = c(0,0,0,0))
plot(t1, show.tip.label = F, lwd = 0.1)
for (i in 1:length(t1$node.label)) {
  node = i + Ntip(t1)
  if (t1$node.label[i] < 100) {
    nodelabels("", node, frame = "none", pch = 16, col = "red")
  }
}
dev.off()

length(which(t1$node.label < 100))
length(which(t1$node.label < 100)) / length(t1$node.label)
