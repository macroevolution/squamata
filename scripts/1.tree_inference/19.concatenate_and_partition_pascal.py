import re
import argparse
import glob
import subprocess
import os
import operator

if 'pascal' in os.getcwd():
	basedir = '/Users/pascaltitle/Dropbox/Oz_Crown_Ages/'
else:
	basedir = '/Users/Sonal/Dropbox/squamate_tree/'

# concatenate and partition

indir = basedir + 'finalAlignments_June2020/renamed/'
seqs = glob.glob("%s*renamed*aln" % indir)
mtc = ['ND1', 'ND2', 'ND4', 'CYTB', 'COI']
mt = ['ND1', 'ND2', 'ND4', 'CYTB', 'COI', 'rRNA_12S', 'rRNA_16S']
framecsv = os.path.join(indir, 'check_frame.csv')
outstem = basedir + 'finalAlignments_June2020/renamed/concatenated'

frames1 = {}
f = open(framecsv, 'r')
for l in f:
	if not re.search('match', l) and not re.search('gap', l):
		d = re.split(',', l.rstrip())
		if d[0] not in frames1:
			frames1[d[0]] = {}
		frames1[d[0]][d[1]] = int(d[2])

frames = {}
for loc, lochash in frames1.items():
	frames[loc] = max(lochash.items(), key=operator.itemgetter(1))[0]

def get_seq(file):
	f = open(file, 'r')
	seq = {}
	id = ''
	for l in f:
		if re.search('>', l.rstrip()):
			id = re.search('>(\S+)', l).group(1)
			id = re.sub('^_R_', '', id)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()
	return seq

allseq = {}
loci = {}

for seqfile in seqs:
	locus = re.sub(basedir + 'finalAlignments_June2020/renamed/', '', seqfile)
	locus = re.sub('_onedirection.*', '', locus)

	seq = get_seq(seqfile)
	# fix frame
	seq2 = {}
	for id, s in seq.items():
		if locus in frames:
			frame = int(frames[locus]) - 1 
		else:
			frame = 0
		seq2[id] = s[frame:]

	for sp, s in seq2.items():
		loci[locus] = len(s)

		if sp not in allseq:
			allseq[sp] = {}
		if locus in allseq[sp]:
		 	print('%s %s %s is a species repped >1' % (locus, id, sp))
		allseq[sp][locus] = s

for locus in loci:
	if loci[locus] % 3 != 0:
		print('%s not in multiple of 3!!!' % locus)

curlen = 0
parts = {	'nDNA, c1': [], 
			'nDNA, c2': [],
			'nDNA, c3': [],
			'mtDNA, c1': [],
			'mtDNA, c2': [],
			'mtDNA, c3': []}
for locus in sorted(loci.keys()):
	# ribosmal!
	if locus in mt and locus not in mtc:
		parts[locus] = ['%s-%s' % (curlen + 1, curlen + loci[locus])] 
	# mito coding
	elif locus in mtc:
		parts['mtDNA, c1'].append('%s-%s\\3' % (curlen + 1, curlen + loci[locus]))
		parts['mtDNA, c2'].append('%s-%s\\3' % (curlen + 2, curlen + loci[locus]))
		parts['mtDNA, c3'].append('%s-%s\\3' % (curlen + 3, curlen + loci[locus]))
	# nuc coding
	else:
		parts['nDNA, c1'].append('%s-%s\\3' % (curlen + 1, curlen + loci[locus]))
		parts['nDNA, c2'].append('%s-%s\\3' % (curlen + 2, curlen + loci[locus]))
		parts['nDNA, c3'].append('%s-%s\\3' % (curlen + 3, curlen + loci[locus]))
	curlen = curlen + loci[locus]

pout = outstem + '.partitions'
p = open(pout, 'w')
for part in parts:
	p.write('%s = %s\n' % (part, ','.join(parts[part])))
p.close()

meta = outstem + '.csv'
out = outstem + '.phy'
o = open(out, 'w')
m = open(meta, 'w')

o.write('%s %s\n' % (len(allseq), sum(loci.values())))
m.write('species,num_loci,tot_seq_length,num_mt_loci\n')
for sp in allseq:
	s = ''
	sploci = allseq[sp]
	splen = sum(loci[x] for x in sploci)
	num_mt = len([loc for loc in sploci if loc in mt])

	m.write('%s,%s,%s,%s\n' % (sp, len(sploci), splen, num_mt))

	for locus in sorted(loci.keys()):
		if locus in allseq[sp]:
			s += allseq[sp][locus]
		else:
			s += '-' * loci[locus]
	o.write('%s %s\n' % (sp, s))
o.close()
m.close()