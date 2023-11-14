
## README for data directory


### squamatesCoreTaxonomy_repDB.csv
Contains all recognized squamate species as per Reptile Database circa 2019, along with links to taxon names in the phylogeny, in the Roll et al. 2017 geographic range dataset, and to the ecological datasets from Meiri et al. 2018. An NA in the columns "tree", "geog" or "eco" means that that reptileDB taxon is not present. An NA in the "repDB" column means that there were a few taxa that we kept in the tree even though there wasn't a clear link to an existing ReptileDB taxon. 

### squamatesCoreTaxonomy_treeTaxa.csv
Similar file to squamatesCoreTaxonomy_repDB.csv but organized more directly around the taxa in the phylogeny. Here, a difference from the other file is that there are no NA values in the tree taxon names. Instead all names are present, but you can see from the "dataType" column that some taxa are imputed. 


### alldat.csv

This is the main dataset, combining taxonomic, diversification and trait data.

__fields (obvious fields ommitted)__

- treename: taxon name (matches phylogeny tip label if in tree)
- inTree: binary 0/1, value of 1 if taxon is in tree and not imputed
- dataType: whether data for taxon comes from genbank, genomic dataset, or if it was placed via imputation
- inTopoConstraint: logical, if TRUE, then taxon is part of the phylogenomic constraint
- mass: mass in grams as per Feldman et al. 2016
- completeSVL: snout-vent length in mm as per Feldman et al. 2016 and expanded via imputation
- centroidLong/centroidLat: centroid coordinates in decimal degrees of geographic range
- minLong/maxLong/minLat/maxLat: latitude/longitude limits of geographic range in decimal degrees
- bio1/bio7/bio12: range-wide averages of bioclim variables 1/7/12
- cmi: range-wide averages of climatic aridity index 
- npp: range-wide averages of net primary productivity
- elev: range-wide averages of elevation in meters
- tri: range-wide averages of terrain roughness index
- rangeSize: size of geographic range in km2
- climPC1-6: scores of a PCA of bioclim 1-19+CMI+NPP
- skullPC1-2: first two PC axes of PCA of multivariate skull shape
- combinedKinesis: reduction of skull kinesis traits into a single ordinal trait
- dietPC1/2: first two PC axes of PCA on log-ratio transformed diet proportions
- geog_xxx: binary 0/1 presence/absence for various geographic regions, based on geographic range polygon intersections
- dr: DR speciation rate based on genetic-only 6885-tip tree
- meanImputedDR: mean DR across 100 trees with imputed taxa
- clads/meanImputedCLADS: same as with DR above
- bamm: BAMM tip speciation rates based on the 6885-tip tree
- climateNicheRate: TR tip rate calculated from climPC1-6
- thermalNicheRate: TR tip rate of bioclim 1
- clade_xxx: binary 0/1 for whether taxon belongs to some clades of interest

**note:** all "AncDist" fields refer to net innovation (ancestral distance) and all trait rates refer to the TR tip rate.