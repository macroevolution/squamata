rm(list = ls())
setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

alldat <- read.csv('./data/alldat.csv')

inno = list(
  # 'chem' = 'chemosensory_index',
  'skull' = 'skullPC1',
  'vert' = 'numberPresacralVert',
  'elong' = 'elongationIndex',
  'diet' = 'dietPC1'
)
xlabs = c("skull PC1", "# presacral vertebrae", "elongation index",
          "diet PC1")
alldat$log_clads = log(alldat$meanImputedCLADS)
figs = vector("list", length(inno))
library(ggplot2)
for (i in 1:length(inno)) {
  
  figs[[i]] = 
    ggplot(alldat, aes_string(inno[[i]], "log_clads")) +
    geom_point(aes(col = as.factor(clade_Serpentes)), alpha = 0.7, pch = 16) +
    xlab(xlabs[i]) + 
    ylab(expression("ln " ~ lambda[CLaDS])) +
    theme_classic() +
    scale_color_manual(values = c('slateblue1', 'coral2')) +
    theme(legend.position = "none")
}

library(cowplot)
abcd = plot_grid(plotlist = figs, ncol = 2, labels = c("A", "B", "C", "D"))
save_plot("manuscript/figures/SupplInfo_FigS27.png", abcd, 
          ncol = 2, base_width = 3, base_height = 4)
