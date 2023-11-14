trait preparation and diversification README
============================================

bestBAMMrun: BAMM analysis with family-level sampling fractions
CLaDS: 
	- primarytree_100_imputations.CLaDS_rates.txt: CLaDS speciation tip rates across 100 trees with additional imputed taxa
	- primarytree_100_imputations.convergence_stats.csv: Gelman statistic and number of generations for each of the 100 CLaDS runs. 
	- pseudoposterior_100_imputations.CLaDS_rates.txt: CLaDS speciation tip rates across 100 trees from the pseudoposterior set.
	- pseudoposterior_100_imputations.convergence_stats.csv: Gelman statistic and number of generations for those 100 CLaDS runs.



dietModeling:
	- phylo: directory with phylogenetic category count modeling code and modeling outputs
	- nonphylo: directory with non-phylogenetic category count modeling code and modeling outputs
	- diet.csv: dietary data (full stomach equivalents) used to infer species-level proportional utilization of prey categories under phylogenetic and non-phylogenetic category count models.
	- diet-matrix.csv: the same data as diet.csv, just in a different format.
	- dietname-to-treename.csv: taxonomic conversion from names in diet data to names used in the phylogeny
	- diet-proportions-phylo.csv: species-level proportional utilization of prey categories inferred under the phylogenetic category count model.
    - diet-proportions-nonphylo.csv: species-level proportional utilization of prey categories inferred under the non-phylogenetic category count model.
    - diet-breadth.R: calculation of diet breadth from estimated diet proportions. Diet breadth (using the phylo proportions) appears in ../alldat.csv.
    - diet-pca.R: principal components analysis of estimated diet proportions. PC1 and PC2 (using the phylo proportions) appear in ../alldat.csv.


2D_Adult_Skull_Coord_treenames.csv: Multivariate Procrustes coordinates for skull shape

climateSpace.rds: R data object for Supplemental fig S21. Contains a list of climate values for 10k random locations across land surfaces, and a list of which of these points fall in each species' geographic range polygon.
