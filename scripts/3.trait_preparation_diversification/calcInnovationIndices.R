# calculate innovation / ancestral distances

library(bm) # remotes::install_github('blueraleigh/bm')

setwd('~/Dropbox/Oz_Crown_Ages')

source('./scripts/empiricalScripts/6.generateFullDataset.R')
source('./scripts/empiricalScripts/calcNetInnovation.R')

outfolder <- './trait-data/netInnovation/'

phy = read.newick("./final-trees/best_ultrametric_fulltree_ddBD_revision.tre")

###################################################
## CLIMATE

x <- alldat[, grep('climPC', colnames(alldat), ignore.case = TRUE, value = TRUE)]
rownames(x) <- alldat$treename
x <- x[complete.cases(x), ]

ancdist <- calcNetInnovation(x, phy)

write.csv(ancdist, paste0(outfolder, "climate-ancdist.csv"))



########################################################
# SKULL SHAPE

x = read.csv("./trait-data/dasilva-nc-2018/raw/2D_Adult_Skull_Coord_treenames.csv")
goodcol = grep('Proc', colnames(x), value = TRUE,ignore.case = TRUE)
rownames(x) <- x$treename
x <- x[, goodcol]

ancdist <- calcNetInnovation(x, phy)

write.csv(ancdist, paste0(outfolder, "skull-shape-ancdist.csv"))


########################################################
# DIET 


x = read.csv('./trait-data/diet-data/data-derived/sysbio-method/diet-proportions-r1.csv', row.names=1)
x = data.matrix(x)
x = sweep(log(x), 1, rowMeans(log(x)))

ancdist <- calcNetInnovation(x, phy)

write.csv(ancdist, paste0(outfolder, "diet-anc-dist.csv"))

    
###########################################
# DIET BREADTH

x = read.csv('./trait-data/diet-data/data-derived/sysbio-method/diet-breadth-r1.csv', row.names=1)
x <- setNames(x[,1], rownames(x))

ancdist <- calcNetInnovation(x, phy)

write.csv(ancdist, paste0(outfolder, "dietBreadthAncDist.csv"))


###########################################
# VERTEBRAL COUNTS

x = read.csv("./trait-data/morph-data/squamates-vertebrae-data.csv")
x <- setNames(x$numberPresacralVert, x$tree)
storage.mode(x) <- "double"

hist(x, breaks='fd') # can probably get away with normal approximation
hist(log(x), breaks='fd')

ancdistLn <- calcNetInnovation(log(x), phy)
ancdist <- calcNetInnovation(x, phy)

write.csv(ancdist, paste0(outfolder, "nPresacralVertAncDist.csv"))
write.csv(ancdistLn, paste0(outfolder, "nPresacralVertAncDist_ln.csv"))


###########################################
# LENGTH AND MASS

x = read.csv("./trait-data/SVLimputation/lengthMassComplete.csv")

ancdistSVL <- calcNetInnovation(log(setNames(x$completeSVL, x$treename)), phy)
ancdistMASS <- calcNetInnovation(log(setNames(x$mass, x$treename)), phy)
ancdistELONG <- calcNetInnovation(log(setNames(x$elongationIndex, x$treename)), phy)

write.csv(ancdistSVL, paste0(outfolder, "ancDist_svl.csv"))
write.csv(ancdistMASS, paste0(outfolder, "ancdist_mass.csv"))
write.csv(ancdistELONG, paste0(outfolder, "ancdist_elong.csv"))


