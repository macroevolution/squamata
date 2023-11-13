import re
import glob
import numpy as np
import os
import subprocess
import argparse
import sys

## three bad seqeunces
## inds: too short (<150, <150, <200, <300)
## inds: too underrepresented (<5%, <10%, <30%, <50%)
## columns with too much missing (<5%, <10%, <30%, <50%)
## mike's 

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--seqlen", required=True, help='minimum seq len to retain.')
parser.add_argument("-n", "--perloci", required=True, help='minimum percent loci for ind to retain.')
parser.add_argument("-m", "--miss", required=True, help='minimum completeness to retain.')
parser.add_argument("-d", "--dir", required=True, help='aligment dir.')
parser.add_argument("-o", "--out", required=True, help='outfile stem')
args = parser.parse_args()

seqlendrop = int(args.seqlen)
underrepresented = float(args.perloci)
missing = float(args.miss)
aln_dir = args.dir
mainout = args.out

outdir = re.sub('alignments', 'trim_%s_%s_%s' % (seqlendrop, underrepresented, missing), aln_dir)
if not os.path.isdir(outdir):
	os.mkdir(outdir)
mainout = mainout + '_%s_%s_%s.out' % (seqlendrop, underrepresented, missing)
mo = open(mainout, 'w')

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

	return seq

alns = glob.glob(aln_dir + '/*aln')

badinds = ['UMMZ_200750',  'UMMZ_200756', 'NA_DLR0623_Va_caud']
badloci = {}
inds = {}

for ix, aln in enumerate(alns):
	locname = re.search('([^/]+)\.fasta', aln).group(1)
	seqaln = get_seq(aln)
	for ind in seqaln:
		if ind not in inds:
			inds[ind] = 0
		inds[ind] += 1
	print(ix)

for bad in badinds:
        mo.write('%s,ALL,USER_DELETED\n' % bad)

# inds: too underrepresented
num_loci = len(alns)
for ind in inds:
	percent = inds[ind] / float(num_loci)
	if percent < underrepresented:
		badinds.append(ind)
		mo.write('%s,ALL,TOO_UNDERREPRESENTED\n' % (ind))

# columns with too much missing
for ix, aln in enumerate(alns):
        locname = re.search('([^/]+)\.fasta', aln).group(1)
        seqaln = get_seq(aln)

	for badind in badinds:
		if badind in seqaln:
			del seqaln[badind]

	for ind, s in seqaln.items():
		slen = len(s) - s.count('-') - s.count('N') - s.count('?')
		if slen < seqlendrop:
			del seqaln[ind]
			mo.write('%s,%s,TOO_SHORT\n' % (ind, locname))

	if len(seqaln) > 3:
		seq = [True] * len(list(seqaln.values())[0])
		# for each column
		for ix, pos in enumerate(seq):
			# identify all sites
			sites = [s[ix] for c, s in seqaln.items()]
			# identify portion missing
			permiss = len(sites) - (sites.count('-') + sites.count('N') + sites.count('?'))
			permiss = permiss / float(len(sites))
			# identify if it is good or not
			if permiss < missing:
				seq[ix] = False
		out = os.path.join(outdir, '%s.fasta.aln' % locname)
		o = open(out, 'w')
		for c, s in seqaln.items():
			newseq = ''.join([bp if keep else '' for bp, keep in zip(s, seq)])
			newseq = newseq.upper()
			o.write('>%s\n%s\n' % (c, newseq))
		o.close()
		mo.write('%s,%s,NUM_SITES_REMOVED\n' % (seq.count(False), locname))

		# do mike's
		out2 = os.path.join(outdir, '%s.fasta2.aln' % locname)
		ambi = subprocess.check_output("python ~/2017_SqCL/ambi_sites_filter_GD_SS.py -v %s %s" % (out, out2), shell=True)
		# num seq changed
		change = re.search('changed\s(\d+)', ambi).group(1)
		change = int(change)
		subprocess.call("mv %s %s" % (out2, out), shell=True)
	
	        if change > 0:
        	        mo.write("%s,%s,AMBIGUOUS_FILTERED\n" % (locname, change))
	else:
		mo.write("%s,ALL,TOO_FEW_INDS\n" % (locname))
