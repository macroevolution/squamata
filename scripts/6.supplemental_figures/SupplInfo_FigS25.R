setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(patchwork)

#########################################################
# no variable explains speciation rate variation
#########################################################
sp = readRDS("pgls.pseudoposterior.Rds")
sp2 = as.data.frame(do.call("rbind", sp))
sp3 = gather(sp2, key = "model", value = "r2")

models = c("M_KINESIS",
           "M_MERISTIC2", "M_PARITY",
           "M_MORPH", "M_TROPHIC1",
           "M_TROPHIC3", "M_MERISTIC1",
           "M_TROPHIC2", "M_BIOCLIM",
           "M_SKULL", "M_GEO")
mnames = c("skull kinesis",
           "vertebral count", "parity mode", 
           "body size/shape",  "diet",  "chemosensory", 
           "digit/limb count",  "foraging mode", 
           "climate",  "skull shape", "geography")
names(mnames) = models

sp3$model2 = mnames[ sp3$model ]
pgls_graph = ggplot(sp3, aes(model2, r2)) + geom_boxplot() + xlab("") +
  ylab(expression("adjusted" ~ R^2)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8)) +
  ylim(0, 0.1)

#########################################################
# innovation and rates are higher in snakes than lizards
#########################################################
d = read.csv("./manuscript/figure_data/tr_innovation.pseudoposterior.csv")
alldat <- read.csv('./data/alldat.csv')
alldat[alldat$clade_Serpentes == 1, "clade_Serpentes"] = "snake"
alldat[alldat$clade_Serpentes == 0, "clade_Serpentes"] = "lizard"
d$clade_Serpentes = alldat[match(d$treename, alldat$treename), "clade_Serpentes"]

diet = read.csv("./manuscript/figure_data/diet_tr_inno.no_phylo.csv")
diet$clade_Serpentes = alldat[match(diet$X, alldat$treename), "clade_Serpentes"]
diet_tr = pull(diet %>% group_by(clade_Serpentes) %>% 
                 select(clade_Serpentes, dietRate) %>% 
  summarize_all(mean) %>% 
  spread(clade_Serpentes, dietRate) %>% mutate(diff = snake / lizard) %>%
  select(diff))

mtr = alldat %>% group_by(clade_Serpentes) %>% 
  select(dietRate, elongationRate, skullRate, massRate) %>%
  summarize_all(mean, na.rm = T) %>% gather(type, value, -clade_Serpentes) %>%
  ungroup() %>%
  spread(clade_Serpentes, value) %>% mutate(diff = snake / lizard) %>%
  mutate(metric = gsub("Rate", "", type))
mtr[mtr$metric == "diet", "diff"] = diet_tr

tr = d %>% select(metric, tr, treenum, clade_Serpentes) %>% 
  group_by(treenum, clade_Serpentes, metric) %>% 
  summarize(tr = mean(tr)) %>% ungroup() %>% 
  spread(clade_Serpentes, tr) %>% mutate(diff = snake / lizard)
tr[tr$metric == "diet_noTree", "metric"] = "diet"
tr_graph = ggplot(tr, aes(diff)) + geom_histogram() + 
  geom_vline(data = mtr, aes(xintercept = diff), col = "red") +
  theme_classic() + facet_wrap(~metric, scales = "free_x", ncol = 4) +
  xlab(expression(TR[snakes] ~ "/" ~ TR[lizards])) 

diet_psi = pull(diet %>% group_by(clade_Serpentes) %>% select(clade_Serpentes, ancDistDiet) %>% 
                 summarize_all(mean) %>% 
                 spread(clade_Serpentes, ancDistDiet) %>% mutate(diff = snake / lizard) %>%
                 select(diff))
mpsi = alldat %>% group_by(clade_Serpentes) %>% 
  select(ancDistDiet, ancDistElongationIndex, skullAncDist, ancDistMass) %>%
  summarize_all(mean, na.rm = T) %>% gather(type, value, -clade_Serpentes) %>%
  ungroup() %>%
  spread(clade_Serpentes, value) %>% mutate(diff = snake / lizard)
mpsi$metric = c("diet", "elongation", "mass", "skull")
mpsi[mpsi$metric == "diet", "diff"] = diet_psi
psi = d %>% select(metric, psi, treenum, clade_Serpentes) %>% 
  group_by(treenum, clade_Serpentes, metric) %>% 
  summarize(psi = mean(psi)) %>% ungroup() %>% 
  spread(clade_Serpentes, psi) %>% mutate(diff = snake / lizard)
psi[psi$metric == "diet_noTree", "metric"] = "diet"

psi_graph = ggplot(psi, aes(diff)) + geom_histogram() + 
  geom_vline(data = mpsi, aes(xintercept = diff), col = "red") +
  theme_classic() + facet_wrap(~metric, scales = "free_x", ncol = 4) +
  xlab(expression(psi[snakes] ~ "/" ~ psi[lizards])) 

#########################################################
# statistical approaches identify shifts in 
# traits & rates in early snake nodes (r2)
#########################################################

r = read.csv("./manuscript/figure_data/r2.pseudoposterior.csv")
r = r[r$trait %in% c("skull", "diet_noTree", "mass",
                     "elongation", "climate"), ]
tips = unique(c(r$tip1, r$tip2))

t = ape::read.tree("./data/1.tree_inference/mainTrees/best_ultrametric_fulltree_ddBD_revision.tre")
tips2 = unique(c(tips, sample(t$tip.label, 100)))
t2 = ape::keep.tip(t, tips2)
traits = unique(r$trait)
edgecols = rep("slateblue1", length(t2$edge.length))

sn = t2$tip.label[ t2$tip.label %in% alldat[alldat$clade_Serpentes == "snake", "treename"] ]
al = t2$tip.label[ t2$tip.label %in% alldat[alldat$clade_Alethinophidia == 1, "treename"] ]
snakeNode = phytools::findMRCA(t2, sn)
alNode = phytools::findMRCA(t2, al)
edgecols[which(t2$edge[,1] %in% c(snakeNode, phangorn::Descendants(t2, snakeNode, type = 'all')))] = 'red'
edgecols[which(t2$edge[,1] %in% c(alNode, phangorn::Descendants(t2, alNode, type = 'all')))] = 'orange'

plot_r2_trees <- function(type = "raw") {
  r1 = r[r$type == type, ]
  
  nodes = vector("list", length(traits))
  
  for (i in 1:length(traits)) {
    r2 = r1[r1$trait == traits[i], ]
    nodes[[i]] = list(rep(0, ape::Nnode(t2)), 
                      rep(0, ape::Nnode(t2)))
    for (j in 1:nrow(r2)) {
      tips = r2[j, c("tip1", "tip2")]
      node = phytools::findMRCA(t2, as.character(tips)) - ape::Ntip(t2)
      nodes[[i]][[1]][node] = nodes[[i]][[1]][node] + 1
      nodes[[i]][[2]][node] = nodes[[i]][[2]][node] + r2[j, "r2"]
    }
    nodes[[i]][[2]] = nodes[[i]][[2]] / nodes[[i]][[1]]
  }
  
  par(mfrow = c(1, length(traits)), mar = c(0, 0, 2, 0))
  for (i in 1:length(nodes)) {
    # color edges by snake and non-snake
    ape::plot.phylo(t2, show.tip.label = F, edge.width = 0.5,
               edge.color = edgecols)
    # opacity by # of times, size by r2
    size = ifelse(is.finite(nodes[[i]][[2]]), nodes[[i]][[2]], 0.01)
    txt = ifelse(nodes[[i]][[1]] == 0, "", nodes[[i]][[1]])
    ape::nodelabels("", frame = "none", pch = 16, 
               col = scales::alpha("black", 0.5),
               cex = size * 5)
    ape::nodelabels(txt, frame = "none")
    # need a legend?
    txt = traits[i]
    if (txt == "diet_noTree") {
      txt = "diet"
    }
    mtext(txt, side = 3, line = 0.3)
  }
}

plot_r2_trees("raw")
# plot_r2_trees("trait")

#########################################################
# speciation rates are higher in snakes than lizards
#########################################################
snakes <- alldat[alldat$clade_Serpentes == "snake", 'treename']
lizards <- alldat[alldat$clade_Serpentes == "lizard", 'treename']

maindiff = alldat %>% group_by(clade_Serpentes) %>% 
  summarize(rate = mean(meanImputedCLADS, na.rm = T)) %>%
  spread(clade_Serpentes, rate) %>% mutate(diff = snake / lizard)

# read in rates from pseudo-posterior
rr = read.csv("data/3.trait_preparation_diversification/CLaDS/pseudoposterior_100_imputations.CLaDS_rates.txt",
              sep = " ", row.names = 1)
rr$treename = rownames(rr)
rr1 = rr %>% gather(tree, rate, -treename)
rr1$clade_Serpentes = alldat[ match(rr1$treename, alldat$treename), "clade_Serpentes"]
rr2 = rr1 %>% group_by(clade_Serpentes, tree) %>% 
  summarize(rate = mean(rate)) %>% spread(clade_Serpentes, rate) %>%
  ungroup() %>% mutate(diff = snake / lizard)

lambda_graph = ggplot(rr2, aes(diff)) + geom_histogram(bins = 50) +
  xlab(expression(lambda[snakes] ~ "/ " ~ lambda[lizards])) +
  theme_classic() + 
  geom_vline(xintercept = maindiff$diff, col = "red")

library(cowplot)
ab = plot_grid(lambda_graph, pgls_graph, ncol = 2, align = "h", 
               labels = c("D", "E"))
allgraph = plot_grid(tr_graph, psi_graph, 
          plot_r2_trees, ab, ncol = 1, 
          labels = c("A", "B", "C", ""))
save_plot("./manuscript/figures/SupplInfo_FigS25.pdf", allgraph, base_width = 8, 
          base_height = 10)
