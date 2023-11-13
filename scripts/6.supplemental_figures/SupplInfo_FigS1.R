rm(list = ls())
library(ape)

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

x = read.csv("./data/6.supplemental_figures/phylogenomic_samples.csv")

########
# get ladderized tree for plotting

t1 = read.tree("./data/0.phylogenomic_data/all.ind0.01_loci0.05_n5185.examl.tre")
t2 = root(t1, "taeGut2", resolve.root = TRUE)
t3 = drop.tip(t2, c("chrPic1", "hg38", "H20145a", 
                    "allMis1", "galGal5", "taeGut2",
                    "UMFS-10956c", "ISIS373002"))
t4 = read.tree(text = write.tree(ladderize(t3)))

########
# get data matrix for plotting
d = read.csv("./data/6.supplemental_figures/SupplInfo_FigS1.sample_matrix.csv", header=T, row.names=1)
m = as.matrix(d)
m = m[t4$tip.label,]

########
# matrix includes loci dropped bc of data quality (high lnl levels) 
# and because of low missingness, need to remove
l = read.csv("./data/6.supplemental_figures/SupplInfo_FigS1.drop_loci_lnl.csv", header = T)
counts = apply(m, 2, sum)
# need to filter out loci sampled for < 0.05 inds
drop = names(counts[counts / dim(m)[1] < 0.05])
drop = unique(c(drop, l$locus))
drop = gsub("\\-", ".", drop)
loci = names(counts[order(counts, decreasing=T)])
loci = loci[!loci %in% drop]

########
# reorder loci by completeness
fac = ifelse(rownames(m) %in% x[x$type2 == "new", "sample"], 2, 1)
m = m * fac
m = t(m)
m = m[loci,]

dim(m)
mean(apply(m, 1, function(x) sum(x > 0)) / ncol(m))

edgecols = rep("black", length(t4$edge.length))
tipnums = which(t4$tip.label %in% x[x$type2 == "new", "sample"])
edgecols[ which(t4$edge[,2] %in% tipnums) ] = "#1f78b4"

pdf("manuscript/figures/SupplInfo_FigS1.pdf", width=6.3, height=4)
par(mfrow=c(1,2))
par(mar = c(3.5, 0, 0, 2))
plot.phylo(t4, direction = "rightwards", show.tip.label = FALSE, 
           lwd = 0.2, y.lim = c(1, length(t4$tip.label)),
           edge.color = edgecols)
par(mar = c(4, 0, 0.5, 0.5))
gl <- 1:(length(t4$tip.label) + 1)
image(1:(dim(m)[1] + 1), 1:(length(t4$tip.label) + 1), m, 
      col = c("lightgray", "black", "#1f78b4"), xlim = c(1, dim(m)[1] - 1), 
      ylim = c(1, length(t4$tip.label) - 1), axes=F, xlab="")
xvals = pretty(c(1,dim(m)))[1:6]
axis(1, at=xvals, labels=xvals)
mtext("number of loci", 1, line=2.5)
dev.off()
