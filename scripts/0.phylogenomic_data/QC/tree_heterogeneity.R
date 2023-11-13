library(ape)
library(phangorn)
library(phytools)

dir = '/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/internal_rtrees_300_0.05_0.3/'

trees = list.files(dir, pattern='SH.tre', full.names=T)
loci = gsub("^.*\\/", '', trees)
loci = gsub(".SH.tre", '', loci)

res = data.frame(locus=loci, mean_resid=rep(NA, length(loci)), sum_resid=rep(NA, length(loci)), stringsAsFactors=F)
for (i in 1:length(trees)) {
	t = read.tree(trees[i])

	root_tip = diag(vcv.phylo(t))
	tips = seq(1, length(t$tip.label))
	nodes = sapply(tips, function(x) {length(Ancestors(t, x, type="all"))})
	names(nodes) = t$tip.label
	r = lm(root_tip ~ nodes[names(root_tip)])
	res[i, "mean_resid"] = mean(abs(residuals(r)))
	res[i, "sum_resid"] = sum(abs(residuals(r)))
}

write.csv(res, "/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/likelihood/residual.csv", row.names=F)
