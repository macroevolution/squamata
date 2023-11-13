import re
import pandas as pd
import os
import numpy as np
import argparse
import glob

parser = argparse.ArgumentParser()
# gene tree dir
parser.add_argument("-g", "--gdir", required=True, help='directory with gene trees')
# stem 1
parser.add_argument("-s1", "--stem1", required=True, default="tree1", help='stem 1 name')
# stem2
parser.add_argument("-s2", "--stem2", required=True, default="tree2", help='stem 2 name')
# alndir
parser.add_argument("-i", "--indir", required=True, help='aln dir')
# out
parser.add_argument("-o", "--out", required=True, help='out file')
args = parser.parse_args()

def get_ll(infile):
	if os.path.isfile(infile):
		f = open(infile, 'r')
		lnl = np.nan
		for l in f:
			if re.search('^Final GAMMA', l):
				lnl = re.search('([0-9|\-|\.]+)$', l).group(1)
	else:
		lnl = np.nan
	return lnl

indir = args.indir
files = glob.glob('%s/*%s*tre' % (indir, args.stem1))

loci = [re.search('^.*/([^/]+).%s' % args.stem1, x).group(1) for x in files]
d1 = pd.DataFrame(data={'locus': loci})
d1['%s_ll' % args.stem1] = [None] * len(loci)
d1['%s_ll' % args.stem2] = [None] * len(loci)
d1['gene_ll'] = [None] * len(loci)
d1['tips'] = [None] * len(loci)
d1['seqlen'] = [None] * len(loci)

def get_seq(x):
        f = open(x, 'r')
        id = ''
        seq = {}
        for l in f:
                if re.search('^>(\S+)', l.rstrip()):
                        id = re.search('>(\S+)', l.rstrip()).group(1)
                        id = re.sub('^_R_', '', id)
                        seq[id] = ''
                else:
                        seq[id] += l.rstrip()
	f.close()

	numseq = len(seq)
	seqlen = len(seq.values()[0])

	return numseq, seqlen

for ix, row in d1.iterrows():
	alndir = re.sub('trees', 'trim', args.gdir)
	aln = os.path.join(alndir, '%s.fasta.aln' % row['locus'])
	numseq, seqlen = get_seq(aln)

	d1.ix[ix, 'tips'] = numseq
	d1.ix[ix, 'seqlen'] = seqlen
	d1.ix[ix, '%s_ll' % args.stem1] = get_ll(os.path.join(indir, 'RAxML_info.%s_%s' % (row['locus'], args.stem1)))
	d1.ix[ix, '%s_ll' % args.stem2] = get_ll(os.path.join(indir, 'RAxML_info.%s_%s' % (row['locus'], args.stem2)))
	g_ll = os.path.join(indir, 'RAxML_info.%s_%s' % (row['locus'], 'gene'))
	d1.ix[ix, 'gene_ll'] = get_ll(g_ll)

outfile = args.out	
d1.to_csv(outfile, index=False)
