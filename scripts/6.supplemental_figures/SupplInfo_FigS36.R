rm(list = ls())
setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

library(ape)
library(ggplot2)

# get rates from sampling fractions
load("data/3.trait_preparation_diversification/CLaDS/CLaDS.sampling_fraction/best_ultrametric_fulltree_ddBD_revision.CLaDS_samplingfrac.Rdata")
frac = CladsOutput$lambdatip_map
names(frac) = CladsOutput$tree$tip.label

# get imputed rates - 100 imputations of the same tree
d = read.table("./data/3.trait_preparation_diversification/CLaDS/primarytree_100_imputations.CLaDS_rates.txt",
             sep = " ")
d$sampling_frac = frac[match(d$V1, names(frac))]

# correlation of imputed to not
a = ggplot(d, aes(V2, sampling_frac)) +
  geom_point(pch = 16, alpha = 0.6, size = 0.5) +
  theme_classic() + 
  xlab(expression(lambda[CLaDS] ~ ", imputed")) +
  ylab(expression(lambda[CLaDS] ~ ", sampling fraction")) +
  scale_x_log10() + scale_y_log10()

# error relative to sampling frac
alldat <- read.csv('./data/alldat.csv')
ff = alldat %>% group_by(family) %>% 
  summarize(sp = n(), imputed = sum(dataType == "imputed")) %>% 
  ungroup() %>% mutate(missing = imputed / sp)
miss = setNames(ff$missing, ff$family)

d$mean_clads = apply(d[2:101], 1, mean)
d$family = alldat[match(d$V1, alldat$treename), "family"]
d$miss = miss[ d$family ]
d$ratio = (d$mean_clads - d$sampling_frac) / d$mean_clads

dd = d %>% group_by(family) %>% 
  summarize(miss = mean(miss), ratio = mean(abs(ratio))) %>% ungroup()
c = ggplot(dd, aes(miss, ratio)) +
  geom_point(pch = 16) + xlab("% of family imputed") +
  ylab(expression(frac("imputed - sampling fraction", imputed))) +
  # ylab(expression(frac((imputed - sampl. frac.) , imputed))) +
  theme_classic() 

# difference across imputations
b = ggplot(d, aes(V2, V4)) +
  geom_point(pch = 16, alpha = 0.6, size = 0.5) +
  theme_classic() + 
  xlab(expression(lambda[CLaDS] ~ ", imputed (sample 1)")) +
  ylab(expression(lambda[CLaDS] ~ ", imputed (sample 2)")) +
  scale_x_log10() + scale_y_log10()

# correlation across imputations
corrs = data.frame(t(combn(100, 2)), cor = NA)
for (x in 1:nrow(corrs)) {
  corrs[x, "cor"] = cor.test(d[ , corrs$X1[x] + 1], d[ , corrs$X2[x] + 1])$estimate
}
mean(corrs$cor)
range(corrs$cor)

# correlation of sampling frac to imputations
corrs2 = rep(NA, 100)
for (i in 1:100) {
  corrs2[i] = cor.test(d[ , (i + 1)], d$sampling_frac)$estimate
}
mean(corrs2)
range(corrs2)

abc = plot_grid(a, b, c, ncol = 3, labels = c("A", "B", "C"))
save_plot("manuscript/figures/SupplInfo_FigS36.png", abc, 
          ncol = 3, base_height = 3, base_width = 3)
