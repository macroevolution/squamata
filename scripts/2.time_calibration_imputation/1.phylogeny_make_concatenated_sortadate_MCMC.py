import numpy as np
import re
import os
import pandas as pd
import glob
import argparse
import random

parser = argparse.ArgumentParser()
# sortadate file
parser.add_argument('-s', "--sortadate", required=False, help="sortadate file")
# alndir
parser.add_argument('-a', "--alndir", required=False, help="directory with alignments")
# tree
parser.add_argument('-t', "--taxa", required=False, help="taxa file")
# outfile
parser.add_argument('-o', "--outdir", required=False, help="output file")
# number of loci 
parser.add_argument("-n", "--number", required=False, help="number of loci desired; leave blank if all")
args = parser.parse_args()

def get_seq(aln):
	f = open(aln, 'r')
	tmpseq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			id = re.sub('^_R_', '', id)
			tmpseq[id] = ''
		else:
			tmpseq[id] += l.rstrip()
	f.close()

	loclen = len(list(tmpseq.values())[0])

	return tmpseq, loclen


def get_tips(taxafile):
	d = pd.read_csv(taxafile)
	tips = d['sample'].dropna().tolist()

	return tips

tips = get_tips(args.taxa)
print(len(tips))
d = pd.read_csv(args.sortadate)
d['locus'] = [re.sub('\..*$', '', x) for x in d.locus]

# tiplen = []
# for ix, row in d1.iterrows():
#	seq, loclen = get_seq(os.path.join(args.alndir, '%s.fasta.aln' % row['locus']))
#	numind = len([ind for ind in seq if ind in tips])
#	tiplen.append(numind) 
# d1['tips'] = tiplen
# d1.to_csv("~/Desktop/tipcount.csv")
d = d.dropna()
d = d[d.tips >= np.percentile(d.tips, [90])[0]]
d = d[d.support >= np.percentile(d.support, [90])[0]]
d = d[d['var'] <= np.percentile(d['var'], [10])[0]]
d = d.sort_values(by = ['length'], ascending=[0])
loci = d.locus.tolist()
print(len(loci))


if len(loci) > int(args.number):
	loci = loci[0:int(args.number)]

# partdir = re.sub('trim', 'partitions', args.alndir)
# parts = {'nc': '', 'p1': '', 'p2': '', 'p3': '', 'uce': ''}
totlen = 0

allseq = {}
for ind in tips:
	allseq[ind] = ''
parts = {}
# go through each aln
for ix, locus in enumerate(loci):
	seq, loclen = get_seq(os.path.join(args.alndir, '%s.fasta.aln' % locus))
	# parts = get_partition(partdir, locus, totlen, parts)

	parts[locus] = [totlen + 1, totlen + loclen]
	totlen += loclen

	for ind in allseq:
		if ind in seq:
			allseq[ind] += seq[ind]
		else:
			allseq[ind] += '-' * loclen

logf = os.path.join(args.outdir, 'concat_n%s.out' % (len(loci)))
lf = open(logf, 'w')
drop_inds = []
for ind in allseq:
	if re.search('^[\-|N|n]+$', allseq[ind]):
		lf.write('dropping %s\n' % ind)
		drop_inds.append(ind)

numinds = len(allseq) - len(drop_inds)
lf.write('number of loci: %s\n' % len(loci))
lf.write('number of inds: %s\n' % numinds)

for ix, i in enumerate(loci):
	lf.write('locus %s: %s\n' % (ix, i))
lf.close()


if not os.path.isdir(args.outdir):
	os.mkdir(args.outdir)
out1 = os.path.join(args.outdir, 'concat_n%s.phy' % (len(loci)))
out2 = os.path.join(args.outdir, 'concat_n%s.partitions' % (len(loci)))

o = open(out1, 'w')
pout = open(out2, 'w')

o.write('%s\t%s\n' % (numinds, totlen))
# make alignment
for ind in allseq:
	if ind not in drop_inds:
		o.write('%s\t%s\n' % (ind, allseq[ind]))
o.close()

for key in parts:
	pout.write('DNA, %s=%s-%s\n' % (key, parts[key][0], parts[key][1]))
	# if len(parts[key]) > 0:
	#	pout.write('DNA, %s=%s\n' % (key, parts[key]))
pout.close()