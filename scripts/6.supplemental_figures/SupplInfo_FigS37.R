library(ape)
library(scales)

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

get_discord_heights <- function(t1, t2) {
  hts = rep(NA, Nnode(t1))
  discord = rep(FALSE, Nnode(t1))
  for (i in 1:length(discord)) {
    nodenum = i + Ntip(t1)
    desc = t1$tip.label[ phangorn::Descendants(t1, nodenum, type = "tips")[[1]] ]
    nodenum2 = phytools::findMRCA(t2, desc)
    hts[i] = phytools::nodeheight(t2, nodenum2)
    desc2 = t2$tip.label[ phangorn::Descendants(t2, nodenum2, type = "tips")[[1]] ]
    diff = c(setdiff(desc, desc2), setdiff(desc2, desc))
    if (length(diff) > 0) {
      discord[i] = TRUE
    }
  }
  return(list(discord, hts))
}

setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")
tt = read.tree("./data/1.tree_inference/pseudoposterior/pseudoposterior.100.trees")
alldat <- read.csv('./data/alldat.csv')

t = read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
tlet = t
tlet$node.label = rep("", Nnode(tlet))
x = read.csv("results/divergenceNodeLabels.csv")
for (i in 1:nrow(x)) {
  tips = as.character(x[i, c("spanningTaxon1", "spanningTaxon2")])
  nodenum = phytools::findMRCA(t, tips) - Ntip(t)
  tlet$node.label[nodenum] = x[i, "symbol"]
}
gen = get_subsample(alldat, t, "subfamily")
subt = keep.tip(t, gen)
subtlet = keep.tip(tlet, gen)

tt2 = vector("list", length(tt))
for (i in 1:length(tt)) {
   print(i)
   subtt = keep.tip(tt[[i]], gen)
   tt2[[i]] = get_discord_heights(subt, subtt)
}

dd = do.call("rbind", lapply(tt2, function(x) x[[1]]))
subt$node.label = 1 - (apply(dd, 2, sum) / 100)

# define nodebars
nb = data.frame(min = rep(NA, Nnode(subt)),
                max = rep(NA, Nnode(subt)),
                median = rep(NA, Nnode(subt)),
                num = seq(1, Nnode(subt)))

for (i in 1:Nnode(subt)) {
  desc = subt$tip.label[ phangorn::Descendants(subt, i + Ntip(subt), type = "tips")[[1]] ]
  nodenum = phytools::findMRCA(subt, desc) - Ntip(subt)
  vals = unlist(lapply(tt2, function(x) x[[2]][nodenum]))
  nb[i, "min"] = min(vals)
  nb[i, "max"] = max(vals)
  nb[i, "median"] = median(vals)
}

maxht = max(phytools::nodeHeights(subt))
pdf("manuscript/figures/SupplInfo_FigS37.pdf", width = 8, height = 12)
par(mar = c(3, 0, 0, 8), xpd = TRUE)
plot(subt, show.tip.label = F, lwd = 0.01, 
      edge.color = "black")
lastPP <- get("last_plot.phylo",envir=.PlotPhyloEnv)
for (j in seq(50, 200, 50)) {
  loc = max(lastPP$x.lim) - j
  lines(x = c(loc, loc), y = c(-3, max(lastPP$y.lim)),
        lty = "dotted")
}
poly = c(seq(10, maxht, 10), round(maxht / 10, 0) * 10 + 10)
for (j in 1:length(poly)) {
  if (j %% 2) {
    loc1 = max(lastPP$x.lim) - (poly[j] - 10)
    loc2 = max(lastPP$x.lim) - poly[j]
    polygon(x = c(loc1, loc1, loc2, loc2), 
            y = c(-3, max(lastPP$y.lim), max(lastPP$y.lim), -3),
            col = alpha("black", 0.05),
            border = alpha("black", 0.05))
  }
}


for (i in 1:nrow(nb)) {
  lines(x = c(nb[i, "min"], nb[i, "max"]),
        y = rep(lastPP$yy[ i + Ntip(subt)], 2),
        lwd = 2, col = alpha("dodgerblue", 0.7))
}

tiplabels(names(gen[match(subt$tip.label, gen)]), cex = 0.75,
          adj = c(0, 0.5), frame = "none")
nodelabels("", frame = "none", pch = 21, bg = "red3",
           cex = ifelse(subt$node.label <= 0.95, 0.7, 0))
nodelabels("", frame = "none", pch = 21, bg = "white",
           cex = ifelse(subtlet$node.label != "", 3, 0))
nodelabels(subtlet$node.label, frame = "none")
axisPhylo()
dev.off()
