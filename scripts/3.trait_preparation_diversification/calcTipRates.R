# Import trait datasets and calculate derived traits, such as uni- or multivariate rates

# Specifically:
## climate niche tip rate
## thermal niche tip rate
## size tip rate
## body elongation tip rate
## limb/digit change tip rate
## vertabral evolution tip rate
## skull shape rate

setwd('~/Dropbox/Oz_Crown_Ages')
f1 <- './trait-data/svlImputation/lengthMassComplete.csv'
f2 <- './trait-data/climateEnv/climateMat.csv'
f3 <- './trait-data/dasilva-nc-2018/raw/2D_Adult_Skull_Coord_treenames.csv'
f4 <- './trait-data/morph-data/squamates-limb-data.csv'
f5 <- './trait-data/morph-data/squamates-vertebrae-data.csv'
f6 <- './trait-data/foragingMode/foragingMode_allsquamates.csv'
f7 <- './trait-data/parityMode_PyronBurbrink2014/squamParity.csv'
f8 <- './trait-data/Schwenk/jawPrehension_chemosensory.csv'
f9 <- './trait-data/morph-data/squamates-kinesis-data.csv'
f10 <- './trait-data/diet-data/data-derived/sysbio-method/diet-proportions-r1.csv'
f11 <- './trait-data/diet-data/data-derived/sysbio-method/diet-breadth-r1.csv'


masterFile <- './taxon-attributes/squamatesCoreTaxonomy_treeTaxa.csv'


treefile <- './final-trees/best_ultrametric_fulltree_ddBD_revision.tre'

outfile <- './trait-data/traitRates.csv'

library(phylo) # remotes::install_github('blueraleigh/phylo')
library(bm) # remotes::install_github('blueraleigh/bm')
library(mk) # remotes::install_github('blueraleigh/mk')

lengthMass 			<- read.csv(f1)
clim 				<- read.csv(f2)
colnames(clim)[1] 	<- 'treename'
skullCoords 		<- read.csv(f3)
limbdat 				<- read.csv(f4)
colnames(limbdat)[colnames(limbdat) == 'tree'] <- 'treename'
vertebraldat 			<- read.csv(f5)
colnames(vertebraldat)[colnames(vertebraldat) == 'tree'] <- 'treename'
foragingMode 		<- read.csv(f6)
parity 				<- read.csv(f7)
chemosensory		<- read.csv(f8)
kinesis				<- read.csv(f9)
colnames(kinesis)[colnames(kinesis) == 'tree'] <- 'treename'
jawPrehension		<- kinesis[, c('treename', 'prehensionMechanism')]

cladesFile <- './trait-data/cladeMemberships.csv'

clades <- read.csv(cladesFile)

# collapse kinesis dataset into categorical variable.
combinedKinesis <- rep(NA, nrow(kinesis))
names(combinedKinesis) <- kinesis$treename
kinesis$mesokinesis <- ifelse(kinesis$mesokinesis == 'yes', 1, 0)
kinesis$metakinesis <- ifelse(kinesis$metakinesis == 'yes', 1, 0)
kinesis$hypokinesis <- ifelse(kinesis$hypokinesis == 'yes', 1, 0)

table(rowSums(kinesis[, c('mesokinesis', 'metakinesis', 'hypokinesis')]))
combinedKinesis[which(rowSums(kinesis[, c('mesokinesis', 'metakinesis', 'hypokinesis')]) == 0)] <- 'akinetic'
combinedKinesis[which(rowSums(kinesis[, c('mesokinesis', 'metakinesis', 'hypokinesis')]) == 1)] <- 'lowKinesis'
combinedKinesis[which(rowSums(kinesis[, c('mesokinesis', 'metakinesis', 'hypokinesis')]) == 2)] <- 'midKinesis'
combinedKinesis[which(rowSums(kinesis[, c('mesokinesis', 'metakinesis', 'hypokinesis')]) == 3)] <- 'highKinesis'

combinedKinesis[clades[which(clades$Alethinophidia == 1), 'treename']] <- 'hyperkinesis'

combinedKinesis <- combinedKinesis[!duplicated(names(combinedKinesis))]

table(combinedKinesis, useNA = 'always')



master <- read.csv(masterFile)

tree <- read.newick(treefile)


# -----------------------
# From length and mass
## rate of body length evolution
## rate of body mass evolution
## rate of body elongation index evolution

commonTaxa <- intersect(lengthMass[!is.na(lengthMass$completeSVL), 'treename'], tiplabels(tree))
traitdat <- setNames(lengthMass[, 'completeSVL'], lengthMass$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
svlRate <- bm.tiprate(log(traitdat), traitTree)
names(svlRate) <- names(traitdat)

plot(svlRate, type = 'h')

commonTaxa <- intersect(lengthMass[!is.na(lengthMass$mass), 'treename'], tiplabels(tree))
traitdat <- setNames(lengthMass[, 'mass'], lengthMass$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
massRate <- bm.tiprate(log(traitdat), traitTree)
names(massRate) <- names(traitdat)

plot(massRate, type = 'h')

commonTaxa <- intersect(lengthMass[!is.na(lengthMass$elongationIndex), 'treename'], tiplabels(tree))
traitdat <- setNames(lengthMass[, 'elongationIndex'], lengthMass$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
elongationRate <- bm.tiprate(log(traitdat), traitTree) # log transform
names(elongationRate) <- names(traitdat)

plot(elongationRate, type = 'h')

# ------------------------
# From climate data
## multivariate climate niche lability
## thermal niche lability


# multivariate climate niche lability based on the 19 bioclim variables + CMI + NPP
# We will use the first 6 PC axes of a PCA, as those explain > 95% of the variance.
clim <- read.csv('./trait-data/climateEnv/climateMat.csv', row.names = 1)
climVar <- 7:27
clim <- clim[, climVar]
clim <- clim[complete.cases(clim),]
climpca <- prcomp(clim, center = TRUE, scale. = TRUE)
summary(climpca)

commonTaxa <- intersect(rownames(climpca$x), tiplabels(tree))
traitTree <- keep.tip(tree, commonTaxa)
# traitdat <- clim[tiplabels(traitTree), grep('bio\\d\\d?', colnames(clim))]
traitdat <- climpca$x[tiplabels(traitTree), 1:6]
rates <- bm.mvtiprate(traitdat, traitTree)

climateNicheRate <- numeric(nrow(traitdat))
names(climateNicheRate) <- rownames(traitdat)
for (ii in 1:length(rates)){
	climateNicheRate[ii] <- sum(diag(rates[[ii]]))
	# eig <- eigen(rates[[ii]])
	# climateNicheRate[ii] <- sqrt(max(eig$values))
}


plot(climateNicheRate, type = 'h')

# thermal niche rate, univariate rate based on mean annual temperature -- bioclim1
traitdat <- clim[commonTaxa, 'bio1']
names(traitdat) <- commonTaxa
traitTree <- keep.tip(tree, commonTaxa)
table(is.na(traitdat))
traitdat <- traitdat[tiplabels(traitTree)]
thermalNicheRate <- bm.tiprate(traitdat, traitTree)
names(thermalNicheRate) <- names(traitdat)

plot(thermalNicheRate, type = 'h')


# for the 3 variables used for overall climate niche rate, calculate the univariate rates as well
# bio12
traitdat <- clim[commonTaxa, 'bio12']
names(traitdat) <- commonTaxa
traitTree <- keep.tip(tree, commonTaxa)
table(is.na(traitdat))
traitdat <- traitdat[tiplabels(traitTree)]
bio12rate <- bm.tiprate(traitdat, traitTree)
names(bio12rate) <- names(traitdat)

plot(bio12rate, type = 'h')

# bio7
traitdat <- clim[commonTaxa, 'bio7']
names(traitdat) <- commonTaxa
traitTree <- keep.tip(tree, commonTaxa)
table(is.na(traitdat))
traitdat <- traitdat[tiplabels(traitTree)]
bio7rate <- bm.tiprate(traitdat, traitTree)
names(bio7rate) <- names(traitdat)

plot(bio7rate, type = 'h')


# NPP
traitdat <- clim[commonTaxa, 'npp']
names(traitdat) <- commonTaxa
traitTree <- keep.tip(tree, commonTaxa)
table(is.na(traitdat))
traitdat <- traitdat[tiplabels(traitTree)]
nppRate <- bm.tiprate(traitdat, traitTree)
names(nppRate) <- names(traitdat)

plot(nppRate, type = 'h')


# -------------------------
# from skull coordinate data
## rate of skull shape evolution

rownames(skullCoords) <- skullCoords$treename
commonTaxa <- intersect(skullCoords$treename, tiplabels(tree))
traitdat <- skullCoords[commonTaxa, grep('ProcCoord', colnames(skullCoords))]
traitTree <- keep.tip(tree, commonTaxa)

goodcol <- paste("ProcCoord", 1:40, sep="")
traitdat <- traitdat[tiplabels(traitTree), goodcol]

rates <- bm.mvtiprate(traitdat, traitTree)

skullRate <- numeric(nrow(traitdat))
names(skullRate) <- rownames(traitdat)
for (ii in 1:length(rates)){
	skullRate[ii] <- sum(diag(rates[[ii]]))
	# eig <- eigen(rates[[ii]])
	# skullRate[ii] <- sqrt(max(eig$values))
}

plot(skullRate, type = 'h')


# -------------------------
## DIET RATE


diet = data.matrix(read.csv(f10, row.names=1))

subtree <- keep.tip(tree, rownames(diet))

dietTR = matrix(0, Ntip(subtree), 4, dimnames=list(
    tiplabels(subtree), c("identity", "log", "alr", "clr")))

alr = log(sweep(diet[,-31], 1, diet[,31], "/"))
clr = sweep(log(diet), 1, rowMeans(log(diet)))

dietTR[, 1] = sapply(bm.mvtiprate(diet, subtree), function(m) sum(diag(m)))
dietTR[, 2] = sapply(bm.mvtiprate(log(diet), subtree), function(m) sum(diag(m)))
dietTR[, 3] = sapply(bm.mvtiprate(alr, subtree), function(m) sum(diag(m)))
dietTR[, 4] = sapply(bm.mvtiprate(clr, subtree), function(m) sum(diag(m)))

# -----------------------------
## DIET BREADTH

dietbreadth <- read.csv(f11)
dietbreadth <- setNames(dietbreadth[,2], dietbreadth[,1])
subtree <- keep.tip(tree, names(dietbreadth))

dietbreadthTR = matrix(0, Ntip(subtree), 2, dimnames=list(
    tiplabels(subtree), c("identity", "log")))

dietbreadthTR[, 1] = bm.tiprate(dietbreadth, subtree)
dietbreadthTR[, 2] = bm.tiprate(log(dietbreadth), subtree)


# --------------------------
# from morphological data
## rate of limb/digit change
## rate of vertebral evolution

# rate of limb change
commonTaxa <- intersect(limbdat[!is.na(limbdat$numberLimbs), 'treename'], tiplabels(tree))
traitdat <- setNames(limbdat[, 'numberLimbs'], limbdat$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
traitdat <- as.factor(traitdat)
table(traitdat)
limbRate <- mk.tiprate(traitdat, traitTree)
names(limbRate) <- names(traitdat)

plot(limbRate, type = 'h')


# rate of digit change
commonTaxa <- intersect(limbdat[!is.na(limbdat$numberDigits), 'treename'], tiplabels(tree))
traitdat <- setNames(limbdat[, 'numberDigits'], limbdat$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
traitdat <- as.factor(traitdat)
table(traitdat)
digitRate <- mk.tiprate(traitdat, traitTree)
names(digitRate) <- names(traitdat)

plot(digitRate, type = 'h')


# rate of vertebral evolution - treat as discrete or continuous?
commonTaxa <- intersect(vertebraldat[!is.na(vertebraldat$numberPresacralVert), 'treename'], tiplabels(tree))
traitdat <- setNames(vertebraldat[, 'numberPresacralVert'], vertebraldat$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
traitdat <- setNames(as.numeric(traitdat), names(traitdat))
table(traitdat)
vertebralRate <- bm.tiprate(traitdat, traitTree)
names(vertebralRate) <- names(traitdat)

plot(vertebralRate, type = 'h')

# --------------------------------
## rate of foraging mode evolution
commonTaxa <- intersect(foragingMode[!is.na(foragingMode$foragingMode), 'treename'], tiplabels(tree))
traitdat <- setNames(foragingMode[, 'foragingMode'], foragingMode$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
traitdat <- as.factor(traitdat)
table(traitdat)
foragingRate <- mk.tiprate(traitdat, traitTree)
names(foragingRate) <- names(traitdat)

plot(foragingRate, type = 'h')

## rate of parity mode evolution
commonTaxa <- intersect(parity[!is.na(parity$combinedParity), 'treename'], tiplabels(tree))
traitdat <- setNames(parity[, 'combinedParity'], parity$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
traitdat <- as.factor(traitdat)
table(traitdat)
parityRate <- mk.tiprate(traitdat, traitTree)
names(parityRate) <- names(traitdat)

plot(parityRate, type = 'h')


## rate of jaw prehension mode evolution
commonTaxa <- intersect(jawPrehension[!is.na(jawPrehension$prehensionMechanism), 'treename'], tiplabels(tree))
traitdat <- setNames(jawPrehension[, 'prehensionMechanism'], jawPrehension$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
traitdat <- as.factor(traitdat)
table(traitdat)
jawPrehensionRate <- mk.tiprate(traitdat, traitTree)
names(jawPrehensionRate) <- names(traitdat)

plot(jawPrehensionRate, type = 'h')


## rate of cranial kinesis evolution
commonTaxa <- intersect(names(na.omit(combinedKinesis)), tiplabels(tree))
traitdat <- combinedKinesis
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
traitdat <- as.factor(traitdat)
table(traitdat)
kinesisRate <- mk.tiprate(traitdat, traitTree)
names(kinesisRate) <- names(traitdat)

plot(kinesisRate, type = 'h')



## rate of chemosensory index evolution -- number of synapomorphies
### we will model this one as increasing only

commonTaxa <- intersect(chemosensory[!is.na(chemosensory$chemosensory_index), 'treename'], tiplabels(tree))
traitdat <- setNames(chemosensory[, 'chemosensory_index'], chemosensory$treename)
traitTree <- keep.tip(tree, commonTaxa)
traitdat <- traitdat[tiplabels(traitTree)]
traitdat <- as.factor(traitdat)
table(traitdat)

cost = matrix(Inf, nlevels(traitdat), nlevels(traitdat))
diag(cost) = 0
for (i in seq_along(levels(traitdat))) {
    xx = as.numeric(levels(traitdat)[i])
    for (j in seq_along(levels(traitdat))) {
        yy = as.numeric(levels(traitdat)[j])
        if (yy > xx)
            cost[i, j] = yy - xx
    }
}

mpr <- mpr.sankoff(traitTree, traitdat, cost)

tiprate = function (mpr, cost, phy) {
    n = Ntip(phy)
    y = numeric(n)
    for (i in 1:n) {
        anc = ancestors(phy)[[i]]
        for (j in 1:length(anc)) {
            s0 = which.min(mpr[anc[j], ])
            s1 = which.min(mpr[children(phy, anc[j])[1],])
            s2 = which.min(mpr[children(phy, anc[j])[2],])
            cont = cost[s0, s1] + cost[s0, s2]
            d = cont/sum(brlens(phy)[children(phy, anc[j])])
            y[i] = y[i] + d/2^j
        }
    }
    y
}

chemosensoryRate <- tiprate(mpr$mpr, cost, traitTree)
names(chemosensoryRate) <- names(traitdat)

plot(chemosensoryRate, type = 'h')


# ---------------------------------
# ---------------------------------

# write to file

ratemat <- data.frame(treename = master$treename)

ratemat <- cbind(ratemat,
					svlRate = 			svlRate[master$treename],
					massRate = 			massRate[master$treename],
					elongationRate = 	elongationRate[master$treename],
					climateNicheRate = 	climateNicheRate[master$treename],
					thermalNicheRate = 	thermalNicheRate[master$treename],
					bio12rate = 		bio12rate[master$treename],
					bio7rate = 			bio7rate[master$treename],
					nppRate = 			nppRate[master$treename],
					skullRate = 		skullRate[master$treename],
					limbRate = 			limbRate[master$treename],
					digitRate = 		digitRate[master$treename],
					vertebralRate = 	vertebralRate[master$treename],
					foragingRate = 		foragingRate[master$treename],
					parityRate = 		parityRate[master$treename],
					jawPrehensionRate = jawPrehensionRate[master$treename],
					kinesisRate = 		kinesisRate[master$treename],
					chemosensoryRate = 	chemosensoryRate[master$treename],
					dietRate = 			dietTR[, 'clr'][master$treename],
					dietBreadthRate = 	dietbreadthTR[, 'log'][master$treename])

head(ratemat)

write.csv(ratemat, file = outfile, row.names = FALSE)

