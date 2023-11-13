rm(list = ls())
setwd('~/Dropbox/Oz_Crown_Ages/dataArchive')

bdir = getwd()

library(ape)
library(cutphylo)

alldat <- read.csv('./data/alldat.csv')
source('./scripts/4.comparative_analyses/model-fits/fxns-find-clades.R')
source("./scripts/4.comparative_analyses/model-fits/fxns-datasets.R")

# location of pseudoposterior of trees
# now found at ./data/2.time_calibration_imputation/pseudoposterior-imputed-trees.zip
treeFiles <- list.files('./taxon-imputation/pseudoposterior-full-trees/', pattern = '\\.tre$', full.names = TRUE)

get_climate_PC <- function(phy) {
  x <- alldat[, grep('climPC', colnames(alldat), ignore.case = TRUE, value = TRUE)]
  rownames(x) <- alldat$treename
  x <- x[complete.cases(x), ]
  
  common = intersect(rownames(x), phy$tip.label)
  
  phy = keep.tip(phy, common)
  
  x = x[rownames(x) %in% phy$tip.label, ]
  
  y = data.matrix(x)
  y = y[phy$tip.label, ]
  return (list(y=y, phy=phy))
}

res = matrix(NA, nrow = 6 * length(treeFiles), ncol = 6)

for (j in 1:length(treeFiles)) {
  
  treenum = gsub(".*xx", "", treeFiles[j])
  treenum = gsub(".tre", "", treenum)
  
  phy = ape::read.tree(treeFiles[j])
  DATASET = list(
    dataset_skull(phy, BASEDIR = bdir),
    dataset_diet_noTree(phy, BASEDIR = bdir),
    dataset_vertebral(phy),
    dataset_mass(phy),
    dataset_elongation(phy),
    get_climate_PC(phy)
  )
  names(DATASET) = c("skull", "diet_noTree", "vertebral",
                     "mass", "elongation", "climate")
  # LOG FALSE, LOG TRUE, LOG TRUE
  # LOG TRUE, LOG TRUE, LOG TRUE
  
  for (i in 1:length(DATASET)) {
    # raw data
    raw = cut.phylo(DATASET[[i]]$phy, DATASET[[i]]$y, 2, method = "arithmetic")
    raw_tip = getSpanningTips(DATASET[[i]]$phy, node = raw$cuts[2])
    raw_r2 = raw$r.squared[2]
    
    # store the data
    res[6 * (j - 1) + i, ] = c(treenum, names(DATASET)[i], "raw", raw_r2, raw_tip)
    cat(i, " * of * ", j, "\n", sep = "")
  }
}

colnames(res) = c("treenum", "trait", "type", "r2", "tip1", "tip2")
write.csv(res, "./manuscript/figure_data/r2.pseudoposterior.csv")
