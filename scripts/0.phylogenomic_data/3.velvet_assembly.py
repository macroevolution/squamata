import pandas as pd
import subprocess
import os

outdir = '/scratch/drabosky_flux/sosi/SqCL_July2017'
file = '/scratch/drabosky_flux/sosi/SqCL_July2017/published_samples.csv'
d = pd.read_csv(file)

name = "assemble"
nodes = 1
cpu = 16
mem = 64
hours = 40

samps = d['sample'].tolist()

for ix, sample in enumerate(samps):
	sh_out = '%s%s.sh' % (name, ix)
	o = open(sh_out, 'w')

	o.write("#PBS -N %s%s\n" % (name, ix))
	o.write("#PBS -M sosi@umich.edu\n")
	o.write("#PBS -A drabosky_flux\n")
	o.write("#PBS -l qos=flux\n")
	o.write("#PBS -q flux\n")
	o.write("#PBS -l nodes=%s:ppn=%s,mem=%sgb\n" % (nodes, cpu, mem))
	o.write("#PBS -l walltime=%s:00:00\n" % hours)
	o.write("#PBS -j oe\n")
	o.write("#PBS -V\n")
	o.write("module load bioperl/1.6.1\n")
	o.write("\n")

	o.write("cd /scratch/drabosky_flux/sosi/SqCL_July2017/velvet_assembly/\n")
	o.write("mkdir $PBS_JOBID\n")
	o.write("cd $PBS_JOBID\n")
	o.write("/home/sosi/bin/VelvetOptimiser/VelvetOptimiser.pl -s 21 -e 91 -x 10 -t 16 -d /scratch/drabosky_flux/sosi/SqCL_July2017/velvet_assembly/%s -f '-shortPaired -fastq.gz -separate /scratch/drabosky_flux/sosi/SqCL_July2017/trim_reads/%s_R1.final.fq.gz /scratch/drabosky_flux/sosi/SqCL_July2017/trim_reads/%s_R2.final.fq.gz -short -fastq.gz /scratch/drabosky_flux/sosi/SqCL_July2017/trim_reads/%s_unpaired.final.fq.gz'\n" % (sample, sample, sample, sample))
	o.write("cp /scratch/drabosky_flux/sosi/SqCL_July2017/velvet_assembly/%s/contigs.fa /scratch/drabosky_flux/sosi/SqCL_July2017/velvet_assembly/%s.fasta\n" % (sample, sample))
	o.close()
	

