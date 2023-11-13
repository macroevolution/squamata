import pandas as pd
import subprocess

# module load sratoolkit/2.8.2-1

d = pd.read_csv("/scratch/drabosky_flux/sosi/SqCL_July2017/published_samples.csv")

for ix, row in d.iterrows():
	run = row['SRR']
	url_template = 'anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/{leading3}/{leading6}/{all}/{all}.sra'
	url = url_template.format(leading3=run[:3], leading6=run[:6], all=run)
	out = '/scratch/drabosky_flux/sosi/SqCL_July2017/prev_published/' + run + '.sra'
	subprocess.call("ascp -k1 -Tr -i /home/sosi/.aspera/cli/etc/asperaweb_id_dsa.openssh %s %s" % (url, out), shell=True)
	out2 = '/scratch/drabosky_flux/sosi/SqCL_July2017/prev_published/' + run + '.fastq'
	# subprocess.call("fastq-dump --split-files %s %s" % (out, out2), shell=True)
