import argparse
import os
import pandas as pd
import re
import subprocess
import zipfile

"""
Sonal Singhal
created on 16 June 2016
Written assuming Trinity 2.2.0
and paired end reads
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Assemble reads using Trinity. "
                            "Written assuming Trinity 2.2.0",
        	formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# trinity program
	parser.add_argument(
                '--trinity',
                type=str,
                default=None,
                help='Full path to Trinity program.'
                )	

	# sample
	parser.add_argument(
                '--sample',
                type=str,
                default=None,
                help='Sample for which to run script.'
                )

	# dir
        parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help="Base directory if running as pipeline."
           	)
           	
        # memory
	parser.add_argument(
                '--mem',
                type=int,
                default=2,
                help='RAM to use for assembly in Gb'
               )

	# CPU
	parser.add_argument(
                '--CPU',
                type=int,
                default=1,
                help='# of CPUs to use in assembly.'
               )

	parser.add_argument(
		'--normal',
		action="store_true",
		default=False,
		help="Run read normalization?."
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

	# read2
        parser.add_argument(
                '--un',
                type=str,
                default=None,
                help="Full path to unpaired file if "
                     "you aren't running in context of pipeline."
                )
                
	# outdir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help="Output directory for Trinity assembly."
		)

	return parser.parse_args()


def get_reads(args):
	# did the reader define the reads?
	if args.read1 == None and args.read2 == None:
		# no she didn't
		read1 = os.path.join(args.dir, 'trim_reads', '%s_R1.final.fq.gz' % args.sample)
		read2 = os.path.join(args.dir, 'trim_reads', '%s_R2.final.fq.gz' % args.sample)
		read_un = os.path.join(args.dir, 'trim_reads', '%s_unpaired.final.fq.gz' % args.sample)
	else:
		read1 = args.read1
		read2 = args.read2
		read_un = args.read_un

	# combine the unpaired and R1 reads as req'd by Trinity
	if read_un != None:
		newread1 = os.path.join(args.outdir, '%s_R1_un.final.gz' % args.sample)
		subprocess.call("cat %s %s > %s" % (read1, read_un, newread1), shell=True)
	else:
		newread1 = read1

	return newread1, read2


def run_trinity(args, read1, read2):
	if not args.outdir:
		outdir = os.path.join(args.dir, 'trinity_assembly')
	else:
		outdir = args.outdir

	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	# trinity requires a dir to exist for it do the assembly
	subdir = os.path.join(outdir, "%s_trinity" % args.sample)

	if not os.path.isdir(subdir):
		os.mkdir(subdir)

	# do this so that there is enough
	# RAM per butterfly
	# cpus = int(args.mem / 10)
	cpus = args.CPU

	if args.normal == False:
		subprocess.call("%s --seqType fq --max_memory %sG --left %s --right %s --CPU %s --output %s" % 
		 		(args.trinity, args.mem, read1, read2, cpus, subdir), shell=True)
	else:
		subprocess.call("%s --normalize_reads --seqType fq --max_memory %sG --left %s --right %s --CPU %s --output %s" %
                                (args.trinity, args.mem, read1, read2, cpus, subdir), shell=True)

	return outdir, subdir


# http://stackoverflow.com/questions/1855095/how-to-create-a-zip-archive-of-a-directory
def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file))


def cleanup(args, read1, outdir, subdir):
	# get rid of read 1 in some cases
	if read1 != args.read1:
		os.remove(read1)

	# make new trinity fasta
	oldfa = os.path.join(subdir, 'Trinity.fasta')
	stem = '%s_contig' % args.sample
	newfa = os.path.join(outdir, '%s.fasta' % args.sample)
	
	id = ''
	seq = {}
	ix = 1	

	f = open(oldfa, 'r')
	o = open(newfa, 'w')

	for l in f:
		if re.search('>', l):
			id = stem + str(ix)
			seq[id] = ''
			ix += 1
		else:
			seq[id] += l.rstrip()
	f.close()

	for id, s in seq.items():
		o.write('>%s\n%s\n' % (id, s))
	o.close()

	# compress outdir
	# dirzip = os.path.join(args.outdir, '%s.zip' % args.sample)
	# zipf = zipfile.ZipFile(dirzip, 'w', zipfile.ZIP_DEFLATED)
	# zipdir(outdir, zipf)
	# zipf.close()


def main():
	# get arguments
	args = get_args()
	# get read names
	# and combine reads if necessary
	read1, read2 = get_reads(args)
	# run trinity
	outdir, subdir = run_trinity(args, read1, read2)
	# cleanup
	cleanup(args, read1, outdir, subdir)
	
if __name__ == "__main__":
	main()
