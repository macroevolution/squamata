import re
import subprocess
import os
import argparse
import pandas as pd

dir = '/scratch/drabosky_flux/sosi/SqCL_July2017/'
mitobim = '/home/sosi/bin/MITObim/MITObim.pl'
runfile = os.path.join(dir, 'squamate_phylogenomics_v3.csv')
mtg = os.path.join(dir, 'metadata', 'reptile_mitogenomes.fa')
mira_dir = '/sw/lsa/centos7/mira/4.0.2/bin/'

parser = argparse.ArgumentParser(description="run for the species",
				 formatter_class=argparse.ArgumentDefaultsHelpFormatter
				)
parser.add_argument('--ind', type=str, default=None, help='individual for which to run.')
parser.add_argument('--outdir', type=str, default=None, help='outdir.')
parser.add_argument('--CPU', type=str, default=None, help='cpu.')
parser.add_argument('--mem', type=str, default=None, help='mem.')
args = parser.parse_args()
ind = args.ind
outdir = args.outdir
cpu = args.CPU
mem = args.mem

d = pd.read_csv(runfile)
out = os.path.join(dir, '%s_mitogenome.fasta' % ind)

def generate_conf(output_dir, individual, ref_path, ref_strain, reads1, reads2, readsu, threads, max_mem):
	output_path = os.path.join(output_dir, "manifest.conf")
	pf_file = open(output_path, 'w')
   	pf_file.write("project = %s\n" % (individual))
	pf_file.write("job = genome,mapping,accurate\n" )
	pf_file.write("parameters = -GE:not=%s -GE:mps=%s -NW:cnfs=warn -NW:mrnl=0 -NW:cac=no -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no\n" % (str(threads), str(max_mem)))
	pf_file.write("\n")
	pf_file.write("readgroup\n")
	pf_file.write("is_reference\n")
	pf_file.write("data = %s\n" % (ref_path))
	pf_file.write("strain = %s\n" % (ref_strain))
	pf_file.write("\n")
	pf_file.write("readgroup = DataIlluminaPairedLib\n")
	pf_file.write("autopairing\n")
	pf_file.write("data = %s %s\n" % (reads1, reads2))
	pf_file.write("technology = solexa\n")
	pf_file.write("template_size = 50 1000 autorefine\n")
	pf_file.write("segment_placement = ---> <---\n")
	pf_file.write("strain = %s\n" % (individual))
	pf_file.write("\n")
	pf_file.write("readgroup = DataIlluminaUnpaired\n")
	pf_file.write("data = %s\n" % (readsu))
	pf_file.write("technology = solexa\n")
	pf_file.write("strain = %s\n" % (individual))
	pf_file.close()
	return output_path

def prep_reads(reads, outdir):
	gzip = []
	orig = []
	for read in reads:
		subprocess.call("gunzip %s" % read, shell=True)
		new = re.sub('.gz', '', read)
		orig.append(new)		
		new1 = re.sub('fq', 'fastq', new)
		subprocess.call("ln -s %s %s" % (new, new1), shell=True)
		gzip.append(new1)
		

	inter_fasta = os.path.join(outdir, '%s.fastq' % ind)
	out = open(inter_fasta, 'w')
	f1 = open(gzip[0], 'r')
	f2 = open(gzip[1], 'r')
	while True:
		a1 = f1.readline().rstrip()
		if not a1: break
		b1 = f2.readline().rstrip()
		a1 = re.search('^(\S+)', a1).group(1)
		b1 = re.search('^(\S+)', b1).group(1)

		if not re.search('/1', a1):
			head1 = re.search('^(\S+)', a1).group(1) + '/1'
			head2 = re.search('^(\S+)', b1).group(1) + '/2'
		else:
			head1 = a1
			head2 = b1

		a2 = f1.readline().rstrip()
		a3 = f1.readline().rstrip()
		a4 = f1.readline().rstrip()

		b2 = f2.readline().rstrip()
		b3 = f2.readline().rstrip()
		b4 = f2.readline().rstrip()

		out.write('%s\n%s\n%s\n%s\n' % (head1, a2, a3, a4))
		out.write('%s\n%s\n%s\n%s\n' % (head2, b2, b3, b4))
	f1.close()
	f2.close()

	f3 = open(gzip[2], 'r')
	for l in f3:
		out.write(l)
	out.close()
	f3.close()

	return inter_fasta, gzip, orig

def run_mira(outdir, conf):
	os.chdir(outdir)
	subprocess.call('%s %s' % ('mira', conf), shell=True)

def run_mitobim(outdir, dir, mitobim, ind, refname, inter_fasta, out_maf, mira_dir):
	os.chdir(outdir)
	subprocess.call('%s -start 1 -end 100 --pair --clean -sample %s -ref %s '
                        '-readpool %s -maf %s --mirapath %s' % (mitobim, ind, refname,
			inter_fasta, out_maf, mira_dir), shell=True)
	d = os.getcwd()
	dirs = [os.path.join(d,o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]
	dirs = [x for x in dirs if re.search('iteration\d+$', x)]
	dirs.sort(key=lambda x: os.path.getmtime(x))
	last = dirs[-1]

	assembly = os.path.join(last, '%s-%s_assembly/%s-%s_d_results/%s-%s_out_%s.unpadded.fasta' % 
		(ind, refname, ind, refname, ind, refname, refname))

	new = os.path.join(dir, 'mt_genomes', '%s_mitogenome.fasta' % ind)
	subprocess.call("cp %s %s" % (assembly, new), shell=True)

	return dirs

# print(d.ix[d['sample'] == ind])
refname = d.ix[d['sample'] == ind, 'mt_genome'].tolist()[0]
ref = os.path.join(outdir, '%s.fa' % refname)
subprocess.call("samtools faidx %s %s > %s" % (mtg, refname, ref), shell=True)

readdir = os.path.join(dir, 'trim_reads')
if not os.path.isdir(outdir):
	os.mkdir(outdir)

reads = [os.path.join(readdir, '%s_R1.final.fq.gz' % ind), 
         os.path.join(readdir, '%s_R2.final.fq.gz' % ind),
	 os.path.join(readdir, '%s_unpaired.final.fq.gz' % ind)]
inter_fasta, gzip, orig_reads = prep_reads(reads, outdir) 
conf = generate_conf(outdir, ind, ref, refname, gzip[0], gzip[1], gzip[2], cpu, mem)
run_mira(outdir, conf)
out_maf = os.path.join(outdir, '%s_assembly' % ind, '%s_d_results' % ind, '%s_out.maf' % ind)
inter_fasta = os.path.join(outdir, '%s.fastq' % ind)
mito_outdirs = run_mitobim(outdir, dir, mitobim, ind, refname, inter_fasta, out_maf, mira_dir)

# clean up! clean up! we all like to clean up!
# remove links, compress read data, remove interleaved
for file in orig_reads:
	subprocess.call("gzip %s" % file, shell=True)
[os.remove(x) for x in gzip]
os.remove(inter_fasta)
# remove reference
# os.remove(ref)
