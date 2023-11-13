rm(list = ls())

library(TreeSim)
library(openxlsx)
library(BAMMtools)
library(geiger)
library(phytools)
library(parallel)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/squamate_tree/")
source("taxon-imputation/missing-tips-fxns.R")

# Change paths
vfulls <- read.tree("phylogenetic_inference/tier2/tier2_100.trees")

cons <- read.xlsx("taxon-imputation/squam-constraints-F-DLR-may21.xlsx")
ug <- unique(cons$genus)

impute_tree <- function(kk) {
  
  vx <- vfulls[[kk]]
  
  
  # Use the overall tree-wide lambda under pure-birth for the messy type 3 constraints:
  lambda <- bd.ms(time = max(branching.times(vx)), n = 10757)
  
  
  for (ii in 1:length(ug)){
    
    tmp <- cons[cons$genus == ug[ii], ]
    cat(ii, "\t||GENUS: ", tmp$genus[1], "\t\t|| ADDING: ", tmp$n_missing[1], "\n", collapse="")	
    
    if (tmp$group[1] == 1 | tmp$group[1] == 2){
      vx <- addMissingTips(vx, tmp$X1[1], tmp$X2[1], 
                           missing=tmp$n_missing[1], in_tree = tmp$n_in_tree[1], 
                           prefix = tmp$genus[1], stem = T)		
    }else if (tmp$group[1] == 3){
      vx <- addMissingTips(vx, tmp$X1[1], tmp$X2[1], missing=tmp$n_missing[1], 
                           in_tree = tmp$n_in_tree[1], stem = T, prefix = tmp$genus[1], lambda=lambda)	
    }else if (tmp$group[1] == 4){
      
      vx <- addMissingTips_wExclusions(vx, cc = tmp, prefix=tmp$genus[1])
      
    }
    
  }
  
  fname <- paste("phylogenetic_inference/tier2/imputations/squams-full-imputed-xx", kk, ".tre", sep="")
  write.tree(vx, file = fname)
}

reps = 1:100
filenames = sapply(reps, function(kk) paste("phylogenetic_inference/tier2/imputations/squams-full-imputed-xx", kk, ".tre", sep=""))
redoreps = reps[!file.exists(filenames)]
mclapply(redoreps, impute_tree, mc.cores = 6)