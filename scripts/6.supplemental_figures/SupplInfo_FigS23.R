rm(list = ls())
library(ape)

setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

alldat <- read.csv('./data/alldat.csv')

t = read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")

fam = unique(alldat$family)
inds = rep(NA, length(fam))
for (i in 1:length(fam)) {
  xx = alldat[alldat$family == fam[i], ]
  xx1 = xx[xx$dataType == "genomic", ]
  xx2 = xx[xx$dataType == "genbank", ]
  if (nrow(xx1) > 0) {
    inds[i] = sample(xx1$treename, 1)
  } else if (nrow(xx2) > 0) {
    inds[i] = sample(xx2$treename, 1)
  }
}
names(inds) = fam
t2 = keep.tip(t, inds)
t2$tip.label = names(inds)[ match(t2$tip.label, inds) ]

xx = alldat %>% group_by(family, dataType) %>% summarize(n = n()) %>% ungroup()
xx1 = xx %>% tidyr::spread(key = dataType, value = n)
xx2 = t(as.matrix(xx1 %>% select(-family)))
colnames(xx2) = xx1$family
xx3 = xx2[, match(t2$tip.label, colnames(xx2))]
xx3[is.na(xx3)] <- 0

xx4 = data.frame(t(xx3))
nrow(xx4[xx4$genomic > 0, ])
nrow(xx4[xx4$genbank > 0, ])
nrow(xx4[xx4$genbank > 0 | xx4$genomic > 0, ])

cols = c("#66c2a5", "#fc8d62", "#8da0cb")
types = c("genomic", "genbank", "imputed")

pdf("manuscript/figures/SupplInfo_FigS23.pdf", width = 8, height = 8)
par(mfrow = c(1, 2))
par(mar = c(4.1, 0, 0, 4.5), xpd = T)
plot.phylo(t2, show.tip.label = F)
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tiplabels(t2$tip.label, cex = 0.7, frame = "none", adj = 0, font = 3)
legend(x = 5, y = 10, fill = cols, legend = types, bty = "n")
par(mar =  c(4.1, 0, 0, 1.1), xpd = T)
plot(NULL, xlim = c(0, max(apply(xx3, 2, sum))),
     ylim = pp$y.lim, axes = F, xlab = "", ylab = "")
for (i in 1:ncol(xx3)) {
  start = 0
  for (j in 1:length(types)) {
    len = xx3[types[j], i]
    end = start + len
    polygon(x = c(start, end, end, start), col = cols[j],
            y = c(i - 0.5, i - 0.5, i + 0.5, i + 0.5), border = NA)
    start = end
  }
}
labels = c(0, 500, 1000, 1500)
axis(1, labels, labels=labels, lwd = 1, line = -0.6, las = 1)
mtext("# of species", side=1, line=1.7)
dev.off()
