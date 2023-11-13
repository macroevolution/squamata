# ASSEMBLE DIVERSIFICATION RATES
## CLADS (with sampling fractions)
## CLADS (fully sampled imputed trees, no sampling fractions)
## BAMM with family-level sampling fractions
## DR from genetic-only tree
## mean DR across fully sampled imputed trees

library(ape)
library(pbapply)

setwd('~/Dropbox/Oz_Crown_Ages')

DR = function(phy) {
    x = numeric(ape::Ntip(phy))
    for (i in 1:ape::Ntip(phy)) {
        p = c(i, ancestors(i))
        for (j in seq_along(p))
            x[i] = x[i] + brlens(p[j]) / 2^(j-1)
    }
    1 / x
}

do_ancestors = function(phy) {
    paths = ape::nodepath(phy)
    function(tip) {
        return (rev(head(paths[[tip]], -1)))
    }
}

do_tips = function(phy) {
    paths = ape::nodepath(phy)
    function(ancestor) {
        which(sapply(paths, function(p) ancestor %in% p))
    }
}

do_brlens = function(phy) {
    len = numeric(ape::Ntip(phy) + phy$Nnode)
    len[phy$edge[,2]] = phy$edge.length
    function(node) {
        len[node]
    }
}

# -------------------------------------------

# Read in phylogeny
tr <- read.tree("./final-trees/best_ultrametric_fulltree_ddBD_revision.tre")
tr <- read.tree(text = write.tree(ladderize(tr)))

# CLADS WITH SAMPLING FRACTIONS
# Load CLADS analysis for tree with sampling fractions (not imputed)
load('./diversification/CLaDS.sampling_fraction/best_ultrametric_fulltree_ddBD_revision.CLaDS_samplingfrac.Rdata')

# make sure tree is the same as the one embedded in the CLADS output
all.equal.phylo(tr, CladsOutput$tree) # TRUE

# get CLADS tip rates
cladstips_samplingFracs <- setNames(CladsOutput$lambdatip_map, tr$tip.label)


# CLADS FROM FULLY SAMPLED IMPUTED TREES
cladsImputed <- data.table::fread('./diversification/primarytree_100_imputations.CLaDS_rates.txt', data.table = FALSE)
rownames(cladsImputed) <- cladsImputed[,1]
cladsImputed <- cladsImputed[, -1]
colnames(cladsImputed) <- NULL
cladsImputed <- as.matrix(cladsImputed)
meancladsImputed <- rowMeans(cladsImputed)

# BAMM WITH SAMPLING FRACTIONS
bammrates <- read.csv('./diversification/bamm/bammTipRates.csv', row.names = 1)
bammrates <- setNames(bammrates[,1], rownames(bammrates))


## DR FROM GENETIC-TAXA-ONLY TREE
brlens = do_brlens(tr)
tips = do_tips(tr)
ancestors = do_ancestors(tr)
tipDR <- DR(tr)
names(tipDR) <- tr$tip.label

plot(tipDR, meancladsImputed[names(tipDR)])
plot(tipDR, bammrates[names(tipDR)])
plot(bammrates, meancladsImputed[names(bammrates)])

## MEAN DR ACROSS FULLY SAMPLED IMPUTED TREES
### set of main imputed trees now found at: data/2.time_calibration_imputation/primarytree-imputed-trees.zip
treeset <- lapply(list.files('./taxon-imputation/best-full-trees', pattern = '\\.tre$', full.names = TRUE), read.tree)
drImputed <- pblapply(treeset, epm::DRstat, cl = 8)
drImputed <- do.call(cbind, drImputed)
meandrImputed <- rowMeans(drImputed)

par(mfrow=c(1,2))
plot(tipDR, meandrImputed[names(tipDR)], xlab = 'tip DR', ylab = 'mean tip DR from imputed', xlim = c(0,2.5), ylim = c(0,2.5))
title(main = 'DR')
abline(a=0, b=1, col = gray(0.5))
plot(cladstips_samplingFracs, meancladsImputed[names(cladstips_samplingFracs)], xlab = 'tip clads', ylab = 'mean tip clads from imputed', xlim = c(0,1.5), ylim = c(0,1.5))
title(main = 'CLADS')
abline(a=0, b=1, col = gray(0.5))

# --------------------------------------------
# Assemble into dataframe
template <- read.csv('./taxon-attributes/squamatesCoreTaxonomy_treeTaxa.csv')
head(template)

divrates <- template[, c('treename', 'family', 'genus', 'dataType')]
divrates <- divrates[order(divrates$family, divrates$genus, divrates$treename), ]
divrates$dr <- tipDR[divrates$treename]
divrates$meanImputedDR <- meandrImputed[divrates$treename]
divrates$bamm <- bammrates[divrates$treename]
divrates$clads <- cladstips_samplingFracs[divrates$treename]
divrates$meanImputedCLADS <- meancladsImputed[divrates$treename]

head(divrates)



write.csv(divrates, './trait-data/diversificationRates.csv', row.names = FALSE)
