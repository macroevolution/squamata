rm(list = ls())
setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

library(ape)
library(ggplot2)
library(phytools)

alldat <- read.csv('./data/alldat.csv')

a = ggplot(alldat, aes(meanImputedCLADS, bamm)) +
  geom_point(pch = 16, alpha= 0.3, aes(col = as.factor(clade_Serpentes))) +
  xlab(expression(lambda[CLaDS])) +
  ylab(expression(lambda[BAMM])) +
  theme_classic() + scale_x_log10() + 
  scale_y_log10() +
  scale_color_manual(values = c('slateblue1', 'coral2')) +
  theme(legend.position = "none")
b = ggplot(alldat, aes(meanImputedCLADS, meanImputedDR)) +
  geom_point(pch = 16, alpha= 0.3, aes(col = as.factor(clade_Serpentes))) +
  xlab(expression(lambda[CLaDS])) +
  ylab(expression(lambda[DR])) +
  theme_classic() + scale_x_log10() + 
  scale_y_log10() +
  scale_color_manual(values = c('slateblue1', 'coral2')) +
  theme(legend.position = "none")
c = ggplot(alldat, aes(bamm, meanImputedDR)) +
  geom_point(pch = 16, alpha= 0.3, aes(col = as.factor(clade_Serpentes))) +
  xlab(expression(lambda[BAMM])) +
  ylab(expression(lambda[DR])) +
  theme_classic() + scale_x_log10() + 
  scale_y_log10() +
  scale_color_manual(values = c('slateblue1', 'coral2')) +
  theme(legend.position = "none")

ab = plot_grid(a, b, c, ncol = 3)
save_plot("manuscript/figures/SupplInfo_FigS43.png", ab, ncol = 3,
          base_height = 3, base_width = 3)
