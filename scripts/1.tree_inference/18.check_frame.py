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

seqs = glob.glob(basedir + "finalAlignments_June2020/renamed/*renamed.aln")
seqs = [seq for seq in seqs if not re.search('rRNA', seq) ]
mt = ['ND1', 'ND2', 'ND4', 'CYTB', 'COI']


def get_seq(seq, file):
	f = open(file, 'r')
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

outfile = 'check_frame.csv'
out = open(outfile, 'w')

for seqfile in seqs:
	locus = re.sub(basedir + 'finalAlignments_June2020/renamed/', '', seqfile)
	locus = re.sub('_onedirection.*', '', locus)
	if re.search('RAG1', locus):
		ref = basedir + 'reference_proteins/RAG1_prot.fa'
	else:	
		ref = os.path.join(basedir + 'reference_proteins/', '%s_prot.fa' % locus)
	print(locus)

	# create one seq hash
	seq = {}
	seq = get_seq(seq, seqfile)
	tmpseq = seqfile + '_2'
	o = open(tmpseq, 'w')
	for id, s in seq.items():
		s = re.sub('-', 'N', s)
		o.write('>%s\n%s\n' % (id, s))
	o.close()

	# double check no frameshifts
	for id, s in seq.items():
		gaps = re.findall('\-+', s)
		for ix, gap in enumerate(gaps):
			if ix == 0 or ix == (len(gaps) - 1):
				pass
			else:
				multi = len(gap) % 3
				if multi != 0:
					out.write('%s,%s,wrong_gap\n' % (locus, id))

	#  do the exonerate
	exonout = ''
	if locus in mt:
		exonout = subprocess.check_output("exonerate -t %s -q %s -m protein2dna "
			"--geneticcode 2 --showalignment F "
			"--showvulgar F --showtargetgff T "
			"| grep 'similarity'" % (tmpseq, ref), shell = True)
	else:
		exonout = subprocess.check_output("exonerate -t %s -q %s -m protein2dna "
			"--showalignment F --showvulgar F --showtargetgff T "
			"| grep 'similarity'" % (tmpseq, ref), shell = True)

	exonout = re.split('\n', exonout)
	starts = {}
	for l in exonout:
		d = re.split('\t', l.rstrip())
		if len(d) > 3:
			start = int(d[3])
			if d[0] not in starts:
				starts[d[0]] = start
			else:
				if start < starts[d[0]]:
					starts[d[0]] = start

	frames = {}
	for id in seq:
		if id in starts:
			frame = starts[id] % 3
			if frame not in frames:
				frames[frame] = 0
			frames[frame] += 1
		else:
			out.write('%s,%s,no_match\n' % (locus, id))
	for frame, ct in frames.items():
		out.write('%s,%s,%s\n' % (locus, frame, ct))

	os.remove(tmpseq)




