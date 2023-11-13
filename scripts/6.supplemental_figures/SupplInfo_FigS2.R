rm(list = ls())
library(ape)
library(ggplot2)
library(tidyverse)
library(cowplot)

setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

d = read.csv("phylogenetic_inference/squamate_phylogenomics_v11_ncbi.csv",
             stringsAsFactors = F)
x = read.csv("./data/6.supplemental_figures/phylogenomic_samples.csv")

# only plot those samples in final tree
x2 = x[which(x$inTopoConstraint == "TRUE"), ]

d1 = d[match(x2$sample, d$sample), ]
d1$tree = x2[match(d1$sample, x2$sample), "type2"]

a = ggplot(d1, aes(avg_length)) + 
  geom_histogram(aes(fill = tree), alpha = 0.3, 
                 position = "identity", bins = 40) +
  xlab("avg. locus length (bp)") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic()
# filter for those with greater
# than zero coverage - individuals that could not be
# assembled new had 0 cov est
b = ggplot(d1 %>% filter(avg_cov > 0), aes(avg_cov)) + 
  geom_histogram(aes(fill = tree), alpha = 0.3, 
                 position = "identity", bins = 40) +
  xlab("avg. coverage") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic()
c = ggplot(d1, aes(num_total)) + 
  geom_histogram(aes(fill = tree), alpha = 0.3, 
                 position = "identity", bins = 40) +
  xlab("# of loci") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic()
# filter for those with greater
# than zero heterozygosity - individuals that could not be
# assembled new had 0 het estimate
d = ggplot(d1 %>% filter(heterozygosity > 0), aes(heterozygosity)) + 
  geom_histogram(aes(fill = tree), alpha = 0.3, 
                 position = "identity", bins = 40) +
  xlab("heterozygosity") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic()

prow <- plot_grid(
  c + theme(legend.position="none"),
  a + theme(legend.position="none"),
  b + theme(legend.position="none"),
  d + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C", "D"),
  hjust = -1,
  nrow = 2
)
prow
legend <- get_legend(
  a + theme(legend.box.margin = margin(0, 0, 0, 12),
            legend.title = element_blank())
)
qual = plot_grid(prow, legend, rel_widths = c(3, .5))
save_plot("manuscript/figures/SupplInfo_FigS2.pdf", qual, 
          base_width = 8, base_height = 5)

# summary stats
d1 %>% filter(tree == "new") %>% 
  select(avg_cov, avg_length, num_total, heterozygosity) %>%
  summarize_all(mean)
table(d1$tree)
