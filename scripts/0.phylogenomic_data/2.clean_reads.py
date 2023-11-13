import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 14 June 2016

Written assuming Trimmomatic v0.36
And PEAR v0.9.10
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Trim Illumina reads for adaptor contamination "
			    "low quality bases and merge reads.\nWritten assuming "
			    "Trimmomatic v0.36 and Pear v0.9.10.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# trimmomatic jar
	parser.add_argument(
                '--trimjar',
                type=str,
            	default=None,
                help='Full path to Trimmomatic jar. Note '
		      'the script assumes you are passing in a'
                      ' .jar file'
                )

	# pear location
	parser.add_argument(
                '--PEAR',
                type=str,
                default=None,
                help='Full path to PEAR binary.'
                )

	# output dir
        parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help='Base directory for pipeline.'
                )

	# sample
	parser.add_argument(
                '--sample',
                type=str,
                default=None,
                help='Sample for which to run script.'
                )

	# sample
	parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with sample info.'
                )

	# CPUs
        parser.add_argument(
                '--CPU',
                type=int,
                default=1,
                help='# of CPUs to use when running Trimmomatic.'
               )

	# min length
	parser.add_argument(
                '--minlength',
                type=int,
                default=36,
                help='Remove reads shorter than this after'
                     ' trimming.'
                )

	# trailing cut
	parser.add_argument(
                '--trail',
                type=int,
                default=3,
                help='Remove trailing low quality bases '
                     'below quality given.'
                )

	# heading cut
	parser.add_argument(
                '--head',
                type=int,
                default=3,
                help='Remove leading low quality bases '
                     'below quality given.'
                )

	# sliding window
	parser.add_argument(
                '--qual',
                type=int,
                default=15,
                help='Cut the read when the 4-base sliding '
		     'window drops below this quality.'
               )
               
        # output dir
	parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Output directory for reads to use '
                     'if not using pipeline.'
                )

	return parser.parse_args()


def sample_info(args):
	# get file data
	d = pd.read_csv(args.file)

	# select data for sample of interest; turn to dict
	d = d.ix[d['sample'] == args.sample]
	d = d[d.read1.notnull()]
	row = d.to_dict('list')

	# convert values from single-item list to values
	info = dict([(key, value[0]) for key, value in row.items()])

	if args.outdir:
		outdir = args.outdir
	else:
		outdir = os.path.join(args.dir, 'trim_reads')

	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	return info, outdir


def rev_comp(seq):
	seq = seq.upper()
	seq_dict = {'A':'T','T':'A','G':'C','C':'G', 'N': 'N'}
	return "".join([seq_dict[base] for base in reversed(seq)])


def get_adapter(a, b, o, aname):
	if re.search('\*', a):
		if pd.isnull(b):
			sys.exit('There is a barcode, but I do not'
			         ' know where to insert it! Add asterisk'
                                 ' to adapter seqeunce where it belongs.')
		else:
			a = re.sub('\*', b, a)
	o.write('>%s\n%s\n>%s_rc\n%s\n' % (aname, a, aname, rev_comp(a)))	
	

def adaptor_file(args, info, dir):
	a_file = os.path.join(dir, '%s_adapters.fa' % args.sample)

	o = open(a_file, 'w')
	get_adapter(info['adaptor1'], info['barcode1'], o, 'ad1')
	get_adapter(info['adaptor2'], info['barcode2'], o, 'ad2')
	o.close()

	return a_file
	

def run_trimmomatic(args, info, a_file, dir):
	out_stem = os.path.join(dir, args.sample)
	
	outfiles = ['%s_R1_paired_1.fq.gz' % out_stem, 
                     '%s_R1_unpaired_1.fq.gz' % out_stem,
                     '%s_R2_paired_1.fq.gz' % out_stem,
                     '%s_R2_unpaired_1.fq.gz' % out_stem]

	subprocess.call("java -jar %s PE -threads %s -phred33 %s %s %s ILLUMINACLIP:%s:2:30:10"  %
                        (args.trimjar, args.CPU, info['read1'], info['read2'], ' '.join(outfiles), 
                        a_file), shell=True)
	 
	return outfiles


def run_pear(args, info, outfiles1, dir):
	out_stem = os.path.join(dir, args.sample)

	pear_out =  ['%s.unassembled.forward.fastq' % out_stem,
		     '%s.unassembled.reverse.fastq' % out_stem,
	             '%s.assembled.fastq' % out_stem,
		     '%s.discarded.fastq' % out_stem]

	outfiles2 = ['%s_R1_paired_2.fq.gz' % out_stem,
                     '%s_R2_paired_2.fq.gz' % out_stem,
                     '%s_R1_assembled_2.fq.gz' % out_stem,
                     '%s_R2_discarded_2.fq.gz' % out_stem]

	subprocess.call("%s -f %s -r %s -o %s -j %s" % (args.PEAR, outfiles1[0], 
                         outfiles1[2], out_stem, args.CPU), shell=True)

	
	for old, new in zip(pear_out, outfiles2):
		subprocess.call("gzip %s" % old, shell=True)
		os.rename(old + '.gz', new)

	return outfiles2


def run_trimmomatic_clean(args, info, outfiles1, outfiles2, dir):
	out_stem = os.path.join(dir, args.sample)

	# combine all single-end read files
	out = '%s_unpaired_2.fq.gz' % out_stem
	single = [outfiles1[1], outfiles1[3], outfiles2[2]]
	subprocess.call("cat %s > %s" % (' '.join(single), out), shell=True)

	outfilesPE = ['%s_R1_paired_3.fq.gz' % out_stem,
                     '%s_R1_unpaired_3.fq.gz' % out_stem,
                     '%s_R2_paired_3.fq.gz' % out_stem,
                     '%s_R2_unpaired_3.fq.gz' % out_stem]	

	# do paired end trimming
	subprocess.call("java -jar %s PE -threads %s -phred33 %s %s %s LEADING:%s "  
                        "TRAILING:%s SLIDINGWINDOW:4:%s MINLEN:%s" %
                        (args.trimjar, args.CPU, outfiles2[0], outfiles2[1], ' '.join(outfilesPE),
                        args.head, args.trail, args.qual, args.minlength), shell=True)
	
	outfileSE = '%s_unpaired_3.fq.gz' % out_stem

	# do single end trimming
	subprocess.call("java -jar %s SE -threads %s -phred33 %s %s LEADING:%s "
                        "TRAILING:%s SLIDINGWINDOW:4:%s MINLEN:%s" %
                        (args.trimjar, args.CPU, out, outfileSE, args.head, args.trail, 
                         args.qual, args.minlength), shell=True)
	os.remove(out)

	return outfilesPE, outfileSE


def clean_up(args, out1, out2, out3, outSE, dir):
	out_stem = os.path.join(dir, args.sample)

	final = ['%s_R1.final.fq.gz' % out_stem,
                 '%s_R2.final.fq.gz' % out_stem,
                 '%s_unpaired.final.fq.gz' % out_stem]

	os.rename(out3[0], final[0])
	os.rename(out3[2], final[1])

	unpaired = [out3[1], out3[3], outSE]
	subprocess.call("cat %s > %s" % (' '.join(unpaired), final[2]), shell=True)

	for file in out1 + out2 + out3 + [outSE]:
		if os.path.isfile(file):
			os.remove(file)


def main():
	# get arguments
	args = get_args()
	# get sample info
	info, dir = sample_info(args)
	# make adaptor file
	a_file = adaptor_file(args, info, dir)
	# run trimmomatic to remove adaptors
	outfiles1 = run_trimmomatic(args, info, a_file, dir)
	# run pear to merge reads
	outfiles2 = run_pear(args, info, outfiles1, dir)
	# run trimmomatic to clean up low quality
	outfiles3, outfileSE = run_trimmomatic_clean(args, info, outfiles1, outfiles2, dir)
	# clean it all up!
	clean_up(args, outfiles1, outfiles2, outfiles3, outfileSE, dir)

if __name__ == "__main__":
	main()
