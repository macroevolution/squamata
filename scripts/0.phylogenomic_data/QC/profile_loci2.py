import re
import pandas as pd
import os
import numpy as np

d1 = pd.read_csv('/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/likelihood/loci_data.csv')
d2 = pd.read_csv('/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/likelihood/residual.csv')
d3 = pd.read_csv('/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/likelihood/sortadate.txt', header=None, names=['locus', 'root_tip_var', 'tree_len'], sep="\t")

d3['locus'] = [re.sub('.SH.tre', '', x) for x in d3['locus']]

d1 = d1.merge(d2, on="locus", how="outer")
d1 = d1.merge(d3, on="locus", how="outer")
d1['astral_ll'] = [None] * d1.shape[0]
d1['concat_ll'] = [None] * d1.shape[0]
d1['gene_ll'] = [None] * d1.shape[0]

def get_ll(infile):
	if os.path.isfile(infile):
		f = open(infile, 'r')
		lnl = ''
		for l in f:
			if re.search('Final GAMMA  likelihood:', l):
				lnl = re.search('Final GAMMA  likelihood: (\S+)', l).group(1)
	else:
		lnl = np.nan
	return lnl

indir = '/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/likelihood'
for ix, row in d1.iterrows():
	d1.ix[ix, 'astral_ll'] = get_ll(os.path.join(indir, 'RAxML_info.%s_astral' % row['locus']))
	d1.ix[ix, 'concat_ll'] = get_ll(os.path.join(indir, 'RAxML_info.%s_concat' % row['locus']))
	d1.ix[ix, 'gene_ll'] = get_ll(os.path.join(indir, 'RAxML_info.%s_gene' % row['locus']))


d1.to_csv('/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/all_likelihood_data.csv', index=False)
