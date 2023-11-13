import re
import gzip
import os
import numpy as np
import argparse

# "heterozygosity"  "num_AHE_contigs" "num_UCE_contigs" "num_total"      
# "avg_cov"         "avg_length"

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample", required=True, help='sample to do.')
args = parser.parse_args()
sample = args.sample

def get_seq(f):
        seq = {}
        id = ''
        
        f = open(f, 'r')
        for l in f:
                if re.search('>', l):
                        id = re.search('>(\S+)', l).group(1)
                        seq[id] = ''
                else:
                        seq[id] += l.rstrip()
        f.close()

        return seq

qual = {'heterozygosity': 0, 'num_AHE_contigs': 0, 
	'num_UCE_contigs': 0, 'num_total': 0,
	'avg_cov': 0, 'avg_length': 0}

prg = '/scratch/drabosky_flux/sosi/SqCL_July2017/PRG/%s.fasta' % sample
vcf = '/scratch/drabosky_flux/sosi/SqCL_July2017/variants/%s.qual_filtered20.cov_filtered2.vcf.gz' % sample

if os.path.isfile(prg):
	seq = get_seq(prg)
	qual['num_total'] = len(seq)
	lens = [len(y) for x, y in seq.items()]
	qual['avg_length'] = np.mean(lens)
	qual['num_AHE_contigs'] = len([x for x in seq if re.search('AHE', x)])
	qual['num_UCE_contigs'] = len([x for x in seq if re.search('uce', x)])

if os.path.isfile(vcf):
	het = {'diff': 0, 'denom': 0}
	cov = 0
	f = gzip.open(vcf, 'r')
	for l in f:
		if not re.search("^#", l):
			het['denom'] += 1
			d = re.split('\t', l.rstrip())
			geno = re.search('^(\S\S\S)', d[9]).group(1)
			geno = re.split('/', geno)
			if geno[0] != geno[1]:
				het['diff'] += 1
			dp = re.search('DP=(\d+)', l).group(1)
			cov += int(dp)
	qual['heterozygosity'] = het['diff'] / float(het['denom'])
	qual['avg_cov'] = cov / float(het['denom'])

out = '/scratch/drabosky_flux/sosi/SqCL_July2017/metadata/quality/%s.csv' % sample
o = open(out, 'w')
vals = sorted(list(qual.keys()))
for val in vals:
	o.write('%s,%s,%s\n' % (sample, val, qual[val]))
o.close()
