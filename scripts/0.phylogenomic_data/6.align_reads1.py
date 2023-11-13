import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 21 June 2016
Written assuming:
	* samtools 1.3.1
	* GATK 3.6
	* picard 2.4.1
	* bwa 0.7.12
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Align reads to lineage, step 1. Written assuming "
                            " samtools 1.3.1, GATK 3.6, picard 2.4.1, and"
                            " bwa 0.7.12",
        	formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# sample
	parser.add_argument(
                '--sample',
                type=str,
                default=None,
                help='Sample for which to run script.'
                )

	# file
	parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with sample info.'
                )

	# basedir
	parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help="Full path to base dir with reads & assemblies."
                )
        
        # bwa
	parser.add_argument(
                '--bwa',
                type=str,
                default=None,
                help='bwa executable, full path.'
                )

	# samtools
        parser.add_argument(
                '--samtools',
                type=str,
                default=None,
                help='samtools executable, full path.'
                )

	# GATK
        parser.add_argument(
                '--gatk',
                type=str,
                default=None,
                help='GATK executable, full path.'
                )

	# picard
        parser.add_argument(
                '--picard',
                type=str,
                default=None,
                help='picard executable, full path.'
                )
	
	# CPUs
	parser.add_argument(
                '--CPU',
                type=int,
                default=1,
                help='# of CPUs to use in alignment.'
               )

	# memory
        parser.add_argument(
                '--mem',
                type=int,
                default=1,
                help='Memory available, as an int, in terms of Gb.'
               )
               
        # outdir
	parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Output directory for alignments, if'
                     ' running out of pipeline.'
                )

	# read1
	parser.add_argument(
                '--read1',
                type=str,
                default=None,
                help="Full path to read 1 file if "
                     "you aren't running in context of pipeline."
                )

	# read2
	parser.add_argument(
		'--read2',
		type=str,
		default=None,
		help="Full path to read 2 file if "
		     "you aren't running in context of pipeline."
		)

	# read u
        parser.add_argument(
                '--un',
                type=str,
                default=None,
                help="Full path to unpaired file if "
                     "you aren't running in context of pipeline."
                )

	# PRG
	parser.add_argument(
                '--prg',
                type=str,
                default=None,
                help="Full path to pseudoref genome if "
                     "you aren't running in context of pipeline."
                )

	return parser.parse_args()


def get_info(args):
	# gets reads
	if not args.read1 and not args.read2:
		read1 = os.path.join(args.dir, 'trim_reads', '%s_R1.final.fq.gz' % args.sample)
		read2 = os.path.join(args.dir, 'trim_reads', '%s_R2.final.fq.gz' % args.sample)
		un = os.path.join(args.dir, 'trim_reads', '%s_unpaired.final.fq.gz' % args.sample)
		reads = [read1, read2, un]
	else:
		reads = [args.read1, args.read2, args.un]

	# gets the lineage associated with the sample
	d = pd.read_csv(args.file)
	lineage = d.ix[d['sample'] == args.sample, 'lineage'].tolist()[0]

	# gets the prg associated with the lineage
	if not args.prg:
		genome = os.path.join(args.dir, 'PRG', '%s.fasta' % args.sample)
	else:
		genome = args.prg

	# defines the outdir
	if args.outdir:
                outdir = args.outdir
        else:
                outdir = os.path.join(args.dir, 'alignments')

	# makes the outdir
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	
	return reads, lineage, genome, outdir


def prepare_seq(args, genome):
	# does all the prep necessary for the PRG
	# if not os.path.isfile(genome + '.bwt'):
	subprocess.call("%s index %s" % (args.bwa, genome), shell=True)
	# if not os.path.isfile(genome + '.fai'):
	subprocess.call("%s faidx %s" % (args.samtools, genome), shell=True)
	out = re.sub('.fasta$', '.dict', genome)
	# if not os.path.isfile(out):
	subprocess.call("java -jar %s CreateSequenceDictionary R=%s O=%s" % 
                                (args.picard, genome, out), shell=True)


def align_seq(args, r, seq, dir):
	out1a = os.path.join(dir, '%s_1.sam' % args.sample)
	out1as = os.path.join(dir, '%s_1.bam' % args.sample)
	out1b = os.path.join(dir, '%s_2.sam' % args.sample)
	out2a = os.path.join(dir, '%s.mateFixed.bam' % args.sample)
	out2b = os.path.join(dir, '%s.bam' % args.sample)
	out3a = os.path.join(dir, '%s_1.mateFixed.sorted.bam' %  args.sample)
	out3b = os.path.join(dir, '%s_2.mateFixed.sorted.bam' % args.sample)
	out4a = os.path.join(dir, '%s_1.rg.mateFixed.sorted.bam' % args.sample)
        out4b = os.path.join(dir, '%s_2.rg.mateFixed.sorted.bam' % args.sample)
	out5a = os.path.join(dir, '%s_1.dup.rg.mateFixed.sorted.bam' % args.sample)
        out5b = os.path.join(dir, '%s_2.dup.rg.mateFixed.sorted.bam' % args.sample)
	out6 = os.path.join(dir, '%s.dup.rg.mateFixed.sorted.bam' % args.sample)
	intervals = os.path.join(dir, '%s.intervals' % args.sample)
	out7 = os.path.join(dir, '%s.realigned.dup.rg.mateFixed.sorted.bam' % args.sample)
	m_1 = os.path.join(dir, '%s_1.metrics' % args.sample)
	m_2 = os.path.join(dir, '%s_2.metrics' % args.sample)

	# need a tmpdir for when sorting BAM files
	tmpdir = os.path.join(dir, args.sample)
	if not os.path.isdir(tmpdir):
		os.mkdir(tmpdir)

	# align
	subprocess.call("%s mem -t %s %s %s %s > %s" % (args.bwa, args.CPU, seq, r[0], r[1], out1a), shell=True)
	subprocess.call("%s mem -t %s %s %s > %s" % (args.bwa, args.CPU, seq, r[2], out1b), shell=True)
	# fixmate
	# note that had used samtools, appears to not work properly
	subprocess.call("%s view -@ %s -b %s > %s" % (args.samtools, args.CPU, out1a, out1as), shell=True)
	subprocess.call("java -jar %s FixMateInformation I=%s O=%s" % (args.picard, out1as, out2a), shell=True)
	subprocess.call("%s view -@ %s -b %s > %s" % (args.samtools, args.CPU,  out1b, out2b), shell=True)
	# sorted
	subprocess.call("%s sort -@ %s -O bam -o %s -T %s %s" % (args.samtools, args.CPU, out3a, tmpdir, out2a), shell=True)
	subprocess.call("%s sort -@ %s -O bam -o %s -T %s %s" % (args.samtools, args.CPU, out3b, tmpdir, out1b), shell=True)
	# readgroup
	subprocess.call("java -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % 
                          (args.picard, out3a, out4a, args.sample, args.sample, args.sample), shell=True)
        subprocess.call("java -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % 
                          (args.picard, out3b, out4b, args.sample, args.sample, args.sample), shell=True)
	# mark read duplicates
	subprocess.call("java -jar %s MarkDuplicates I=%s O=%s ASO=coordinate METRICS_FILE=%s" % 
                       (args.picard, out4a, out5a, m_1), shell=True)
	subprocess.call("java -jar %s MarkDuplicates I=%s O=%s ASO=coordinate METRICS_FILE=%s" % 
                       (args.picard, out4b, out5b, m_2), shell=True)
	subprocess.call("%s merge -@ %s %s %s %s" % (args.samtools, args.CPU, out6, out5a, out5b), shell=True)
	# indel target
	subprocess.call("%s index %s" % (args.samtools, out6), shell=True)
	subprocess.call("java -Xmx%sg -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s" % 
                        (args.mem, args.gatk, seq, out6, intervals, args.CPU), shell=True)
	# indel realigner
	subprocess.call("java -Xmx%sg -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % 
			(args.mem, args.gatk, seq, out6, intervals, out7), shell=True)
	

	# remove the files
	[os.remove(x) for x in [out1a, out1as, out1b, out2a, out2b, out3a, 
                                out3b, out4a, out4b, out5a, out5b, out6, 
			        out6 + '.bai', intervals]]
	
	# remove the dir
	os.rmdir(tmpdir)


def main():
	# get arguments
	args = get_args()
	reads, lineage, genome, dir = get_info(args)
	# prep sequence
	prepare_seq(args, genome)
	# do the alignments
	align_seq(args, reads, genome, dir)

if __name__ == "__main__":
	main()
