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
vfull <- read.tree("final-trees/best_ultrametric_fulltree_ddBD_revision.tre")

cons <- read.xlsx("taxon-imputation/squam-constraints-F-DLR-may21.xlsx")

# Use the overall tree-wide lambda under pure-birth for the messy type 3 constraints:
lambda <- bd.ms(time = max(branching.times(vfull)), n = 10757)

ug <- unique(cons$genus)

impute_tree <- function(kk) {
  
  vx <- vfull
  
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
  
  fname <- paste("taxon-imputation/best-full-trees/squams-full-imputed-xx", kk, ".tre", sep="")
  write.tree(vx, file = fname)
}

REPS <- 100
reps = seq(1, REPS)

filenames = sapply(reps, function(kk) paste("taxon-imputation/best-full-trees/squams-full-imputed-xx", kk, ".tre", sep=""))
redoreps = reps[!file.exists(filenames)]
mclapply(redoreps, impute_tree, mc.cores = 6)









