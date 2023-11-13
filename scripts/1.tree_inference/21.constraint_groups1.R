library(ape)
library(phytools)
library(phangorn)

if (grepl('pascal', getwd())) {
	setwd('~/Dropbox/Oz_Crown_Ages')
} else {
	setwd('~/Dropbox/squamate_tree')
}

# DNA inds
inds = read.csv("finalAlignments_March2020/renamed/concatenated.csv", stringsAsFactors = F)

# taxonomy
d <- read.csv('fasta_March2020/metaTable2.csv', stringsAsFactors = FALSE)

# constraint tree
tr = read.tree("data/all_bootstraps/all.concat_ind0.01_loci0.05_all_n5185.constraint.tre")
tr$node.label <- NULL
  
# rename tips in genomic tree to update
outgroups <- c('taeGut2', 'ISIS373002', 'UMFS-10956c', 'allMis1', 'chrPic1', 'galGal5', 'hg38', 'H20145a')
outgroups %in% tr$tip.label
outgroups <- intersect(outgroups, tr$tip.label)
tr <- root(tr, outgroup = outgroups)
tr <- ladderize(tr)
tr <- drop.tip(tr, outgroups)

tax <- read.csv('~/Dropbox/Oz_Crown_Ages/sqCL_taxa_repDBmatching.csv', stringsAsFactors=FALSE)
tax$repDBtaxon <- gsub('^\\s+|\\s+$', '', tax$repDBtaxon)
tax$manualMatch <- gsub('^\\s+|\\s+$', '', tax$manualMatch)
# manual matching supercedes repDB auto matching, unless no manual match found
tax[which(tax$manualMatch == '' | is.na(tax$manualMatch)), 'manualMatch'] <- tax[which(tax$manualMatch == '' | is.na(tax$manualMatch)), 'repDBtaxon']
head(tax)
table(tax$manualMatch == '' | is.na(tax$manualMatch))
tax[which(tax$manualMatch == '' | is.na(tax$manualMatch)), 'manualMatch'] <- tax[which(tax$manualMatch == '' | is.na(tax$manualMatch)), 'species']

table(tr$tip.label %in% tax$sample)
tr$tip.label[!tr$tip.label %in% tax$sample]

newLabels <- sapply(tr$tip.label, function(x) tax[which(tax$sample == x)[1], 'manualMatch'])
table(newLabels == '' | is.na(newLabels))

# are there any duplicate taxa?
dups <- newLabels[duplicated(newLabels)]

# there are, but which to drop is unimportant -- the topology is identical.
tr <- drop.tip(tr, names(dups))
newLabels <- newLabels[!names(newLabels) %in% names(dups)]

names(newLabels) <- NULL
newLabels <- gsub('\\s+', '_', newLabels)

tr$tip.label <- newLabels

tr
 
 

# check all species in taxonomy
inds$genus = gsub("_.*$", "", inds$species)
inds$family = d[match(inds$genus, d$genus), "family"]
  
# need to combine some families that are not well represented
fams = table(inds$family)

gekko = c("Eublepharidae", "Diplodactylidae", "Carphodactylidae", "Pygopodidae", "Gekkonidae", "Sphaerodactylidae", "Phyllodactylidae")

scinco = c("Scincidae", "Gerrhosauridae", "Cordylidae", "Xantusiidae")

laterata = c("Bipedidae", "Amphisbaenidae", "Gymnophthalmidae", "Teiidae", "Lacertidae", "Rhineuridae", "Alopoglossidae", "Blanidae", "Cadeidae", "Trogonophidae")

iguania = c("Chamaeleonidae", "Agamidae", "Phrynosomatidae", "Iguanidae", "Dactyloidae", "Leiosauridae", "Polychrotidae", "Hoplocercidae", "Leiocephalidae", "Liolaemidae", "Opluridae", "Corytophanidae", "Crotaphytidae", "Tropiduridae")

anguids = c("Helodermatidae", "Diploglossidae", "Anniellidae", "Anguidae", "Xenosauridae", "Shinisauridae", "Varanidae", "Lanthanotidae")

snakes = c("Acrochordidae", "Aniliidae", "Anomalepididae", "Boidae", "Bolyeridae", "Colubridae", "Cylindrophiidae", "Elapidae", "Gerrhopilidae", "Homalopsidae", "Lamprophiidae", "Leptotyphlopidae", "Loxocemidae", "Pythonidae", "Tropidophiidae", "Typhlopidae", "Uropeltidae", "Viperidae", "Xenopeltidae", "Anomochilidae", "Pareidae", "Psammophiidae", "Xenodermidae", "Xenophidiidae", "Xenotyphlopidae")

splitfams = list(gekko, scinco, laterata, iguania, anguids, snakes)
names(splitfams) = c("gekko", "scinco", "laterata", "iguania", "anguids", "snakes")
                       
  
# confirm all names are in original b/c typos
unlist(splitfams) %in% names(fams)
  
# figure out which fams are missing
fams[ !names(fams) %in% unlist(splitfams)]
  
# need to see if I was over-zealous in grouping
# these aren't in tree but i still grouped in
tfams = unique(inds[match(tr$tip.label, inds$species), "family"])
# okay, i only did that for snakes - that's fine
unlist(splitfams)[!unlist(splitfams) %in% tfams]
# and did I miss any families in my tree?
tfams[!tfams %in% unlist(splitfams)]
  
# check counts in each grouping
lapply(splitfams, function(x) sum(fams[x]))
  
mono = vector("list", length(splitfams))
edgecols = rep("gray", nrow(tr$edge))
cols = c("red", "blue", "green", "orange", "black", "purple")
# need to confirm all groupings are monophyletic
for (i in 1:length(splitfams)) {
	splitfam = splitfams[[i]]
    keepinds =  inds[inds$family %in% splitfam, "species"]
    inds[which(inds$family %in% splitfam), "group" ] = names(splitfams)[i]
    tinds = keepinds[keepinds %in% tr$tip.label]
    anc = getMRCA(tr, tinds)
    tinds2 = tr$tip.label[Descendants(tr, anc, type = "tips")[[1]]]
    
    edgecols[tr$edge[, 2] %in% Descendants(tr, anc, type = "all")] = cols[i]
    
    tfams2 = unique(inds[match(tinds2, inds$species), "family"])
    mono[[i]] = setdiff(tfams2, splitfam)
    # include Sphenodon outgroup in constraint tree
    tt = drop.tip(tr, setdiff(tr$tip.label, c(tinds, 'Sphenodon_punctatus')))
    # unroot the constraint tree
    tt <- unroot(tt)
    # make small constraint trees
    write.tree(tt, paste0('tier1/clades/', names(splitfams)[i], ".tre"))
}
  
# write out groups
write.csv(inds, "tier1/clades/inds_to_groups.csv", row.names = F)
  
# plot groupings
pdf("tier1/clades/monophyly_test.pdf")
plot.phylo(tr, edge.color = edgecols, show.tip.label = F)
dev.off()