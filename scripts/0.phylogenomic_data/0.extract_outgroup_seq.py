import re
import subprocess as sp
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genome", required=True, help='Genome file.')
args = parser.parse_args()

genome = args.genome
out = re.sub('.fa', '_blat.out', genome)

seqs = {}
o = open(out, 'r')
for l in o:
	if not re.search('^#', l):
		d = re.split('\s+', l.rstrip())
		if d[0] not in seqs and float(d[10]) <= 1e-20:
			dir = '+'
			if int(d[8]) > int(d[9]):
				dir = '-'

			if re.search('uce', d[0]):
				start = int(d[8]) - 500
				end = int(d[9]) + 500
			else:
				sub = int(d[6]) - 1
				add = 1500 - int(d[7])
	
				if dir == '+':
					start = int(d[8]) - sub
					end = int(d[9]) + add
				else:
					end = int(d[8]) + sub
					start = int(d[9]) - add
		
			if start < 1:
				start = 1
			
			seqs[d[0]] = {'contig': d[1], 'start': start, 'end': end}
o.close()

faidx = genome + '.fai'
if not os.path.isfile(faidx):
	sp.call("~/bin/samtools-1.3.1/samtools faidx %s" % genome, shell = True)

out = re.sub('.fa', '.fasta', genome)
for s in seqs:
	call = "~/bin/samtools-1.3.1/samtools faidx %s %s:%s-%s >> %s" % (genome, seqs[s]['contig'], seqs[s]['start'], seqs[s]['end'], out)
	sp.call(call, shell=True)
