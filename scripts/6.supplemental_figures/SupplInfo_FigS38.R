rm(list = ls())
setwd("~/Dropbox/Oz_Crown_Ages/dataArchive")

library(readxl)
library(ggplot2)
library(cowplot)

d = read_xlsx("./scripts/6.supplemental_figures/data_for_figures/SupplInfo_Fig_S38.previous_phylogenetic_studies.xlsx")
d$year = as.numeric(gsub(".*\\s", "", d$paper))
d$study = FALSE
d[d$paper == "current study 2023", "study"] = TRUE

a = ggplot(d, aes(genomic_loci, genomic_species)) +
  geom_jitter(pch = 21, aes(bg = study),  size = 2.5) + 
  xlab("# phylogenomic loci") +
  ylab("# phylogenomic species") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gray80", "red")) 
b = ggplot(d, aes(genetic_loci, genetic_species)) +
  geom_jitter(pch = 21, aes(bg = study), size = 2.5) + 
  xlab("# genetic loci") +
  ylab("# genetic species") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gray80", "red"))
c = ggplot(d, aes(genomic_species, genetic_species)) +
  geom_jitter(pch = 21, aes(bg = study), size = 2.5) + 
  xlab("# phylogenomic species") +
  ylab("# genetic species") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gray80", "red"))
ab = plot_grid(a, b, c, ncol = 3, labels = c("A", "B", "C"))
save_plot("manuscript/figures/SupplInfo_FigS37.png", ab,
          ncol = 3, base_height = 3, base_width = 3)
