import re
import glob
import os
import pandas as pd
import subprocess

genes = {'ND4': 'c', 'COI': 'c', 'CYTB': 'c', 'ND1': 'c', 'ND2': 'c', 'rRNA_12S': 'nc', 'rRNA_16S': 'nc'}
seqfiles = glob.glob('/nfs/turbo/lsa-rabosky/Lab/SqCL_July2017/mt_genomes/*fasta')
outdir = '/nfs/turbo/lsa-rabosky/Lab/SqCL_July2017/mt_genomes_anno/'

mdir = '/scratch/drabosky_flux/sosi/phylogeny/'

sfile = os.path.join(mdir, 'squamate_phylogenomics_v11.csv')
d = pd.read_csv(sfile)
sp = {}
for ind, species in zip(d['sample'], d['lineage']):
        sp[ind] = species

'''
for seqfile in seqfiles:
	for gene, gtype in genes.items():
		out = seqfile + '.%s' % gene
		out = re.sub('mt_genomes', 'mt_genomes_anno', out)
		genefile = '/scratch/drabosky_flux/sosi/phylogeny/for_pascal/mtDNA_loci/mitogenomes_%s.fa' % gene
		if not os.path.isfile(out):
			if gtype == 'c':
				subprocess.call("exonerate --model coding2coding --geneticcode 2 %s %s --showalignment no --showtargetgff no --showquerygff yes > %s" % (seqfile, genefile, out), shell=True)
			else:
				subprocess.call("exonerate --model affine:local --geneticcode 2 %s %s --showalignment no --showtargetgff no --showquerygff yes > %s" % (seqfile, genefile, out), shell=True)
'''

for gene in genes:
	outfiles = glob.glob('%s*%s' % (outdir, gene))
	matches = {}

	print(gene)
	for ix, outfile in enumerate(outfiles):
		if ix % 10 == 0:
			print(ix)
		sample = re.sub(outdir, '', outfile)
		sample = re.sub('_mitogenome..*', '', sample)

		seqfile = re.sub('.fasta.*$', '.fasta', outfile)
		seqfile = re.sub('mt_genomes_anno', 'mt_genomes', seqfile)

		species = sp[sample]

		o = open(outfile, 'r')
		for l in o:
			if re.search('similarity', l):
				d = re.split('\t', l.rstrip())
				if species in matches:
					length = matches[species]['end'] - matches[species]['start']
					new_len_diff = (int(d[4]) - int(d[3])) / float(length)
					if new_len_diff >= 0.9 and int(d[5]) > matches[species]['score']:
						matches[species] = {'start': int(d[3]), 'end': int(d[4]), 'seqid': d[0],
									'score': int(d[5]), 'seq': seqfile, 'ind': sample}
				else:
					matches[species] = {'start': int(d[3]), 'end': int(d[4]), 'seqid': d[0],
								'score': int(d[5]), 'seq': seqfile, 'ind': sample}
		o.close()

	out = os.path.join(outdir, '%s.fasta' % gene)
	o = open(out, 'w')
	for species in matches:
		call = subprocess.Popen("samtools faidx %s %s:%s-%s" % (matches[species]['seq'], 
			matches[species]['seqid'], matches[species]['start'], matches[species]['end']), 
			shell=True, stdout=subprocess.PIPE)
		lines = [line.rstrip() for line in call.stdout]
		o.write('>%s@%s %s\n' % (species, matches[species]['ind'], matches[species]['score']))
		for l in lines[1:]:
			o.write(l + '\n')
	o.close()
