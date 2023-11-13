import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 23 June 2016
Written assuming nothing!
"""

def get_args():
	parser = argparse.ArgumentParser(
			description="This creates the files that then get " 
                                    "aligned in the next script.",
           		formatter_class=argparse.ArgumentDefaultsHelpFormatter
			)

	# file
	parser.add_argument(
		'--file',
		type=str,
		default=None,
		help='File with information for phylogeny making.'
		)

	# dir
	parser.add_argument(
		'--dir',
		type=str,
		default=None,
		help='Base directory when used in context of '
                     'pipeline.'
	 	)

	# output dir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help='Output directory for alignments if not '
	             'running in context of pipeline.'
		)

	# all or not
	parser.add_argument(
		"--all",
		action="store_true",
		default=False,
		help="all individuals or just internal?"
		)

	# family level
        parser.add_argument(
                "--family",
                action="store_true",
                default=False,
                help="family level?"
                )

	return parser.parse_args()


def get_files(args):
	# make all the directory structure
	if args.outdir:
		outdir = args.outdir
	else:
		outdir = os.path.join(args.dir, 'phylogeny')

	# main result folder
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	x = pd.read_csv(args.file)
	if args.all:
		d = x.ix[x.keep == True]
	elif args.family:
		d = x.ix[x.family_level == True]
	else:
		ours = ['genome', 'pyron_sqcl2', 'rablab_brazil', 'rablab_rapid', 'rablab_sqcl2']
		d = x.ix[x['type'].isin(ours)]
		d = d.ix[d.keep == True]
	
	genomes = {}
	for l in d['sample'].unique().tolist():
		g1 = os.path.join(args.dir, 'PRG', '%s.fasta' % l)
		g2 = os.path.join(args.dir, 'variable_PRG', '%s.fasta' % l)
		if os.path.isfile(g2):
			genomes[l] = g2
		elif os.path.isfile(g1):
			genomes[l] = g1 
		else:
			print(l)

	return outdir, genomes
	

def get_seq(genomes):
	seqs = {}
	ids = {}

	for lin, file in genomes.items():
		seqs[lin] = {}
		seq = os.path.join(file)
                s = open(seq, 'r')
                id = ''

                for l in s:
                        if re.search('>', l):
                                id = re.search('>(\S+)', l.rstrip()).group(1)
                                seqs[lin][id] = ''
				if id not in ids:
					ids[id] = 0
				ids[id] += 1
                        else:
                                seqs[lin][id] += l.rstrip()
                s.close()

	return seqs, ids


def print_loci(outdir, seq, loci):
	sps = sorted(seq.keys())
	n_sp = len(sps)

	d = os.path.join(outdir, 'locus_data.csv')
	d = open(d, 'w')
	d.write('locus,n_lineages,missingness,length,PICs\n')

	for locus in loci:
		# count how many sps have the locus
		count = sum([1 for sp in sps if locus in seq[sp]])

		# only print out the locus if in 4 sp
		if count >= 4:
			out = os.path.join(outdir, '%s.fasta' % locus)
			o = open(out, 'w')

			for sp in sps:
				if locus in seq[sp]:
					o.write('>%s\n%s\n' % (sp, seq[sp][locus]))
			o.close()

		d.write('%s,%s,%.3f,NA,NA\n' % (locus, count, count / float(n_sp)))

	d.close()


def main():
	# get arguments
	args = get_args()
	# get genome files, make dirs
	outdir, gen = get_files(args)
	# get sequences
	seq, loci = get_seq(gen)
	# print the loci
	print_loci(outdir, seq, loci)	

if __name__ == "__main__":
	main()
