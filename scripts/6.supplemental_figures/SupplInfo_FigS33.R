library(ape)

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

x = read.tree("./data/6.supplemental_figures/SupplInfo_FigS33.phyparts.tre")
t1 = x[[1]]
t2 = x[[2]]
t1$node.label = as.numeric(t1$node.label)
t2$node.label = as.numeric(t2$node.label)

c = data.frame(agree = t1$node.label, conflict = t2$node.label)
c$total =  c$agree + c$conflict
c[!complete.cases(c$total), "total"] = 0

pdf("SupplInfo_FigS33.pdf", height=12, width=8)
par(lwd=0.5, mar=c(0, 2, 0, 5), xpd = T)
plot(t1, show.tip.label =  F)
# cols = c("#91cf60", "#d73027")
cols = c("#8babf1", "#c44601")
tiplabels(t1$tip.label, cex = 0.7, frame = "none", adj = 0, font = 3)
for (i in 1:nrow(c)) {
  if (c[i, "total"] > 0) {
    nodelabels(text = "", node = i + Ntip(t1), adj = c(0, 0.5),
               pie = c[i, c("agree", "conflict")], cex = 0.7, piecol=cols)
  }
}

legend(x=0.02, y = 12,
       legend=c("support", "conflict"), 
       fill = cols, bty="n", cex=1.2)
dev.off()
