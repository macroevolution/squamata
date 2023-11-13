import re
import numpy as np
import os 
import glob
import subprocess as sp
import p4
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

def get_seq(aln):
	seq = {}
	id = ''
	f = open(aln, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			id = re.sub('^_R_', '', id)
			seq[id] = ''
		else:
			seq[id] += l.rstrip().upper()
	f.close()

	return seq

def get_miss(seq, length):
	vals = [s.count('-') + s.count('N') for s in seq.values()]
	miss = np.mean(vals) / float(length)
	return miss

def get_het(seq, length):
	hets = [len(re.findall('[r|y|s|w|k|m|R|Y|S|W|K|M]', s)) for s in seq.values()]
	return np.mean(hets) / float(length)

def get_gc(seq, length):
        vals = [s.count('G') + s.count('C') for s in seq.values()]
        gc = np.mean(vals) / float(length)
        return gc

def get_pic(seq, length):
	pic = 0
	redund = 0

	for i in range(0, length):
		bases = [s[i] for s in seq.values()]
		bases = [bp for bp in bases if bp in ['A', 'T', 'C', 'G']]
		bases = list(set(bases))
		if len(bases) > 1:
			pic += 1
		if len(bases) > 2:
			redund += 1

	return [pic / float(length), redund / float(length)]

def get_tree_support(treefile):
	ro.r('t1 = read.tree("%s")' % treefile)
	return ro.r('mean(as.numeric(t1$node.label), na.rm=T)')[0]

def print_reduced(outdir, seq, locus, treefile):
	# get tips
	tips = [tip for tip in seq]
	ro.r.assign('tips', tips)
	ro.r('tips = unlist(tips)')
	
	# get union
	ro.r('keep = intersect(ctree$tip.label, atree$tip.label)')
	if os.path.isfile(treefile):
                ro.r('gtree = read.tree("%s")' % treefile)
                ro.r('gtree$tip.label = gsub("^_R_", "", gtree$tip.label)')
		ro.r('keep = intersect(keep, gtree$tip.label)')

	# drop extra tips in atree
	ro.r('atree1 = drop.tip(atree, setdiff(atree$tip.label, keep))')
	ro.r('atree1$node.label = NULL')
	ro.r('atree1$edge.length = NULL')
	atree1 = os.path.join(outdir, '%s.astral_species.tre' % locus)
	ro.r('write.tree(atree1, "%s")' % atree1)

	# drop extra trips in cocnat tree
	ro.r('ctree1 = drop.tip(ctree, setdiff(ctree$tip.label, keep))')
	ctree1 = os.path.join(outdir, '%s.concat_species.tre' % locus)
	ro.r('write.tree(ctree1, "%s")' % ctree1)

	# do gene tree
	if os.path.isfile(treefile):
		ro.r('gtree1 = drop.tip(gtree, setdiff(gtree$tip.label, keep))')
		gtree1 = os.path.join(outdir, '%s.gene.tre' % locus)
		ro.r('write.tree(gtree1, "%s")' % gtree1)

	keep = ro.r('keep')
	alnout = os.path.join(outdir, '%s.fasta.aln' % locus)
	o = open(alnout, 'w')
	for ind in keep:
		o.write('>%s\n%s\n' % (ind, seq[ind]))
	o.close()
	
	return atree1, ctree1, alnout
	
importr('ape')
alndir = '/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/internal_trim_300_0.05_0.3/'
alns = glob.glob(alndir + '*aln')
treedir = re.sub('trim', 'trees', alndir)
concat = '/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/internal_concat/ExaML_result.concat_ind0.05_loci0.95_all_n3077_2'
astral = '/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/internal_astral/trees_0.05_0.05_10_100_5e-05_SH.all.tre'
outdir = '/scratch/drabosky_flux/sosi/SqCL_July2017/phylogeny/likelihood/'

ro.r('ctree = read.tree("%s")' % concat)
ro.r('atree = read.tree("%s")' % astral)

res = {}
for aln in alns:
	locus = re.sub('^.*/', '', aln)
	locus = re.sub('.fasta.aln', '', locus)
	res[locus] = {}

	# get seq
	seq = get_seq(aln)

	# get length
	res[locus]['length'] = len(seq.values()[0])
	
	# get completeness
	res[locus]['occupancy'] = len(seq)

	# get missing
	res[locus]['missing'] = get_miss(seq, res[locus]['length'])

	# get heterozygosity
	res[locus]['heterozygosity'] = get_het(seq, res[locus]['length'])

	# GC content
	res[locus]['GC'] = get_gc(seq, res[locus]['length'])

	# pics / saturation
	vals = get_pic(seq, res[locus]['length'])
	res[locus]['PICs'] = vals[0]
	res[locus]['saturation'] = vals[1]

	# get tree
	# get average SH
	treefile = os.path.join(treedir, '%s.SH.tre' % locus)
	if os.path.isfile(treefile):
		res[locus]['SH'] = get_tree_support(treefile)
	else:
		res[locus]['SH'] = 'NA'

	# get average boot
	treefile = os.path.join(treedir, 'RAxML_bipartitions.%s' % locus)
	if os.path.isfile(treefile):
		res[locus]['boot'] = get_tree_support(treefile)
	else:
		res[locus]['boot'] = 'NA'

	# print out reduced alignment
	# print out reduced concat tree
	# print out reduced astral tree
	(atree1, ctree1, aln) = print_reduced(outdir, seq, locus, treefile)	


loc1 = res.keys()[0]
types = sorted([type1 for type1 in res[loc1]])
out = open(os.path.join(outdir, 'loci_data.csv'), 'w')
out.write('locus,%s\n' % (','.join(types)))
for locus in res:
	vals = [str(res[locus][type1]) for type1 in types] 
	out.write('%s,%s\n' % (locus, ','.join(vals)))
out.close()
