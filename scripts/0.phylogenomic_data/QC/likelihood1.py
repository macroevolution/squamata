import re
import numpy as np
import random
import os 
import glob
import subprocess as sp
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import argparse

parser = argparse.ArgumentParser()
# tree 1
parser.add_argument("-t1", "--tree1", required=True, help='tree 1 to evaluate')
# tree 2
parser.add_argument("-t2", "--tree2", required=True, help='tree 2 to evaluate')
# stem 1
parser.add_argument("-s1", "--stem1", required=True, default="tree1", help='stem 1 name')
# stem2
parser.add_argument("-s2", "--stem2", required=True, default="tree2", help='stem 2 name')
# alndir
parser.add_argument("-a", "--alndir", required=True, help='aln dir')
# outdir
parser.add_argument("-o", "--outdir", required=True, help='out dir')
args = parser.parse_args()

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

	seq2 = {}
	for id, s in seq.items():
		id = re.sub('_[1|2]$', '', id)
		seq2[id] = s

	return seq2


def print_reduced(outdir, seq, locus, treefile):
	# get tips
	tips = [tip for tip in seq]
	ro.r.assign('tips', tips)
	ro.r('tips = unlist(tips)')
	
	# get union
	ro.r('keep = intersect(tree1$tip.label, tree2$tip.label)')
	ro.r('keep = intersect(keep, tips)')

	# do gene tree
	if os.path.isfile(treefile):
		ro.r('g = read.tree("%s")' % treefile)
		ro.r('g = drop.tip(g, setdiff(g$tip.label, keep))')
		gtree = os.path.join(outdir, '%s.gene.tre' % (locus))
		class1 = ro.r('class(g)')[0]
        	if class1 == "phylo":
        	        ro.r('write.tree(g, "%s")' % gtree)
        	else:
        	        gtree = ''
	else:
		gtree = ''

	# drop extra tips in atree
	ro.r('tree1a = drop.tip(tree1, setdiff(tree1$tip.label, keep))')
	tree1 = os.path.join(outdir, '%s.%s.tre' % (locus, args.stem1))
	class1 = ro.r('class(tree1a)')[0]
	if class1 == "phylo":
		ro.r('write.tree(tree1a, "%s")' % tree1)
	else:
		tree1 = ''

	# drop extra trips in cocnat tree
	ro.r('tree2a = drop.tip(tree2, setdiff(tree2$tip.label, keep))')
	tree2 = os.path.join(outdir, '%s.%s.tre' % (locus, args.stem2))
        class2 = ro.r('class(tree2a)')[0]
        if class2 == "phylo":
		ro.r('write.tree(tree2a, "%s")' % tree2)
	else:
		tree2 = ''

	len_keep = int(ro.r('length(keep)')[0])
	keep = ro.r('keep')
	alnout = os.path.join(outdir, '%s.fasta.aln' % locus)
	if len_keep > 0:
		o = open(alnout, 'w')
		for ind in keep:
			o.write('>%s\n%s\n' % (ind, seq[ind]))
		o.close()
	
	return gtree, tree1, tree2, alnout
	
importr('ape')
alndir = args.alndir
alns = glob.glob(alndir + '/*aln')
treedir = re.sub('trim', 'trees', alndir)
# tree 1 
ro.r('tree1 = read.tree("%s")' % args.tree1)
# tree 2
ro.r('tree2 = read.tree("%s")' % args.tree2)

outdir = args.outdir
if not os.path.isdir(outdir):
	os.mkdir(outdir)

out = open(os.path.join(outdir, 'commands.txt'), 'w')

print(len(alns))
for aln in alns:
	locus = re.sub('^.*/', '', aln)
	locus = re.sub('.fasta.aln', '', locus)

	treefile = os.path.join(treedir, '%s.SH.tre' % locus)

	# get seq
	seq = get_seq(aln)

	# print out reduced alignment
	# print out reduced tree 1
	# print out reduced tree 2
	(gtree, tree1, tree2, aln) = print_reduced(outdir, seq, locus, treefile)
	if (len(tree1) > 0) and (len(tree2) > 0) and (len(gtree) > 0):
		call1 = 'raxmlHPC-PTHREADS -T 2 -m GTRGAMMA -n %s_%s -s %s -g %s.%s.tre -p %s' % (locus, args.stem1, aln, locus, args.stem1, random.randint(1, 1000))	
		call2 = 'raxmlHPC-PTHREADS -T 2 -m GTRGAMMA -n %s_%s -s %s -g %s.%s.tre -p %s' % (locus, args.stem2, aln, locus, args.stem2, random.randint(1, 1000))
		call3 = 'raxmlHPC-PTHREADS -T 2 -m GTRGAMMA -n %s_gene -s %s -g %s.gene.tre -p %s' % (locus, aln, locus, random.randint(1, 1000))
		out.write(call1 + '\n')
		out.write(call2 + '\n')
		out.write(call3 + '\n')
out.close()
