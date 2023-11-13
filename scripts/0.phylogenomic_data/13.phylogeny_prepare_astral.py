import re
import os
import pandas as pd
import glob
import argparse
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indcutoff", required=True, help='percent of loci needed for ind to keep ind')
parser.add_argument("-l", "--loccutoff", required=True, help='percent of inds needed for loci to keep loci')
parser.add_argument("-t", "--tol", required=True, help='branch length to collapse, recommend 5e-5')
parser.add_argument("-c", "--collapse", required=True, help='collapse nodes with <SH or <boot less than this')
parser.add_argument("-b", "--maxbl", required=True, help='remove terminal branches with br lens longer than maxbl * mean(blen)')
parser.add_argument("--boot", action="store_true", default=False, help='use bootstrap trees instead of SH trees')
parser.add_argument("-a", "--alndir", required=True, help='directory with alignments')
parser.add_argument("-o", "--outdir", required=True, help="output directory for astral trees")
parser.add_argument("-d1", "--dropind", required=False, help="csv file with taxa to drop; first column has names")
parser.add_argument("-d2", "--droploci", required=False, help="csv file with loci to drop; first column has names")
parser.add_argument("--AHE", action="store_true", default=False, help='add flag if you only want AHE loci')
parser.add_argument("--uce", action="store_true", default=False, help='add flag if you only want uce loci')

args = parser.parse_args()
indcutoff = float(args.indcutoff)
loccutoff = float(args.loccutoff)
tol = float(args.tol)
collapse = int(args.collapse)
maxbl = int(args.maxbl)
alndir = args.alndir
outdir = args.outdir

# filter inds with a lot of missing data
def inds_drop(inds, numloci, cutoff):
	to_drop = []
	for ind in inds:
		amt = inds[ind] / float(numloci)
		if amt < cutoff:
			to_drop.append(ind)
	return to_drop

# filter loci that have a lot of missing data
def loci_drop(loci, numinds, cutoff):
        to_drop = []
        for locus in loci:
                amt = loci[locus] / float(numinds)
                if amt < cutoff:
                        to_drop.append(locus)
        return to_drop

def manipulate_gene_tree(tree, drop_inds, tol, collapse, maxbl):
	# drop long branches
	ro.r('outs = c("chrPic1", "allMis1", "galGal5", "taeGut2" , "hg38")')
        ro.r('t1 = read.tree("%s")' % tree)
	ro.r('t1$tip.label = gsub("^_R_", "", t1$tip.label)')
	ro.r('blen = t1$edge.length[1:Ntip(t1)]')
	ro.r('weird = t1$edge[which(blen > mean(blen) * %s), 2]' % maxbl)
	ro.r('weird = t1$tip.label[weird]')
	ro.r('weird = as.vector(na.omit(weird))')
	ro.r('weird = weird[!weird %in% outs]')
	ro.r('t2 = drop.tip(t1, weird)')

	# polytomize any weak nodes
	# print(ro.r('Nnode(t2)'))
        ro.r('t3 = di2multi(t2, %s)' % tol)
	# print(ro.r('Nnode(t3)'))

        # collapse any low support branches
	ro.r('t3$node.label = as.numeric(t3$node.label)')
        # print(ro.r('Nnode(t3)'))
	ro.r('t4 = pruneTree(t3, %s)' % collapse)
	# print(ro.r('Nnode(t4)'))

	# drop inds
	ro.r.assign('dropInds', drop_inds)
	ro.r('dropInds = unlist(dropInds)')
	ro.r('dropInds = dropInds[dropInds %in% t4$tip.label]')
	ro.r('t5 = drop.tip(t4, dropInds)')

	# write tree
	tclass = ro.r('class(t5)')[0]
	if tclass == "phylo":
		t = ro.r('write.tree(t5)')[0]
	else:
		t = ''
	return t

def get_inds(aln):
	inds = []
	f = open(aln, 'r')
	for l in f:
		if re.search('>', l):
			ind = re.search('>(\S+)', l).group(1)
			ind = re.sub('_R_', '', ind)
			inds.append(ind)
	f.close()
	return inds

# first need to create a hash of individuals and loci
def get_data(alns):
	inds = {}
	loci = {}
	alns2 = {}
	
	for ix, aln in enumerate(alns):
		locus = re.search('([^/]+).fasta', aln).group(1)
		alns2[locus] = aln

		tmpinds = get_inds(aln)
		for ind in tmpinds:
			if ind not in inds:
				inds[ind] = 0
			inds[ind] += 1
 							
		loci[locus] = len(tmpinds)

	return alns2, inds, loci

treedir = re.sub('trim', 'trees', alndir)

seqtype = 'all'
if args.AHE:
        seqtype = 'AHE'
        alns = glob.glob(alndir + '/AHE*aln')
elif args.uce:
        seqtype = 'uce'
        alns = glob.glob(alndir + '/uce*aln')
else:
        alns = glob.glob(alndir + '/*aln')

if args.dropind:
        xx = pd.read_csv(args.dropind)
        drop_inds1 = xx.ix[:, 0].tolist()
else:
	drop_inds1 = []
if args.droploci:
	xx = pd.read_csv(args.droploci)
        drop_loci1 = xx.ix[:, 0].tolist()

trees, inds, loci = get_data(alns)
drop_inds2 = inds_drop(inds, len(loci), indcutoff)
drop_inds = list(set(drop_inds1 + drop_inds2))
drop_loci2 = loci_drop(loci, len(inds), loccutoff)
drop_loci = list(set(drop_loci1 + drop_loci2))

kept = [x for x in inds if x not in drop_inds] 
locikept = [x for x in loci if x not in drop_loci]

print('number of inds kept: %s' % (len(kept)))
print('number of loci kept: %s' % (locikept))

if args.boot:
	out = os.path.join(outdir, 'trees_%s_%s_%s_%s_%s_boot.%s.trees' % (indcutoff, loccutoff, collapse, maxbl, tol, seqtype))
else:
	out = os.path.join(outdir, 'trees_%s_%s_%s_%s_%s_SH.%s.trees' % (indcutoff, loccutoff, collapse, maxbl, tol, seqtype))
o = open(out, 'w')

importr('phangorn')
importr('ape')
# go through each aln
# dropping any that are in drop loci
for locus in trees:
	if locus not in drop_loci:
		if not args.boot:
			treefile = os.path.join(treedir, '%s.SH.tre' % locus)
		else:
			treefile = os.path.join(treedir, 'RAxML_bipartitions.%s' % locus)
		if os.path.isfile(treefile): 
			tree = manipulate_gene_tree(treefile, drop_inds, tol, collapse, maxbl)
			if (len(tree) > 0):
				o.write(tree + '\n')
			else:
				print('%s missing_tree\n' % locus)
		else:
			print('%s missing_tree\n' % locus)
o.close()
