
DROPBOX = "/Users/drabosky/Dropbox (University of Michigan)/manuscriptwork/SQCL/bigphylo"
if (grepl('pascal|mgrundler', getwd())) {
	DROPBOX = "~/Dropbox/"
}
BASEDIR = file.path(DROPBOX, "Oz_Crown_Ages")

source("fxns-datasets.R")
source("fxns-likelihood.R")
source("fxns-find-clades.R")
source("fxns-shift-search.R")

TREEFILE = file.path(BASEDIR, "./final-trees/best_ultrametric_fulltree_ddBD_revision.tre")
fulltree = ape::read.tree(TREEFILE)

DATASET = c(
      "diet"
    , "diet_noTree"  
    , "skull"
    , "mass"
    , "vertebral"
    , "elongation"
    , "climate"
)
 

#-------------------------#

# Load subfamily skeleton tree, used
#   to constrain possible shift locations:

sfam_tree <- ape::read.tree("subfamily-pairs-squams-r1.tre")
clademat  <- buildSpanningSet_DF(sfam_tree)


