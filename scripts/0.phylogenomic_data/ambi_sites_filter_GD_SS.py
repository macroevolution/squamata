#!/usr/bin/env python2
##useful

from __future__ import print_function
import re
from Bio import SeqIO,Seq,SeqRecord
from collections import Counter
ambi=re.compile('[MRWSKY][^MRWSKY]*[MRWSKY]([^MRWSKY]*[MRWSKY])+')
nongap=re.compile('[^-?]')
acgt20=re.compile('([^ACGTMRWSKY]*[ACGT]){20,}')



def main():
	import argparse
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"-v","--verbose",default=False,
			action="store_true",
			help="""Output diagnostics"""
		)
	parser.add_argument(
			"infile",
			type=str,
			help="""The input file"""
		)
	parser.add_argument(
			"outfile",
			type=str,
			help="""The output file"""
		)
	args=parser.parse_args()
	worker=LocalFilter(args.verbose)
	write_to_file(worker.get_seqs(args.infile),args.outfile)


class LocalFilter(object):
	def __init__(self,verbose=False):
		from sites_local import SiteCounter
		self.proc=SiteCounter()
		self._v=verbose
		self.delta=0
		self.changed=0
	def get_seqs(self, file_path,file_type="fasta"):
		from sites_local import read_sequences_alt
		self.proc.scan_file(file_path,file_type)
		if self._v: 
			print("Sites Scanned")
			import os.path
			self.file_name=os.path.basename(file_path)
		self.delta=0
		self.changed=0
		return self.trim_all(read_sequences_alt(file_path,file_type))
	def trim_seq(self,tmp_rec):
		src_seq=tmp_rec.aligned_seq_data
		foo=ambi.search(str(src_seq))
		if foo is not None:
#             print(tmp_rec.aln_uid, "Removing ambiguous")
			seq=Seq.Seq(sub_ambi(foo),src_seq.alphabet)
		else:
			seq=src_seq
		if self.proc.has_private_alleles(tmp_rec.aln_uid,4):
#             print(tmp_rec.aln_uid, "Scan Private")
			seq,drop=sub_priv(src_seq,
								 self.proc.get_private_alleles(tmp_rec.aln_uid),seq)#,
#                        tmp_rec.aligned_seq_data.alphabet)
			if len(drop) >0:
				if min(drop)>1 and get_acgt_count(src_seq[0:min(drop)-1])<10:
					seq=subst_N_sub(seq,endpos=min(drop)-1,mutable=True,subst='-')
				if max(drop)<len(seq) and get_acgt_count(src_seq[max(drop):])<10:
					seq=subst_N_sub(seq,pos=max(drop),mutable=True,subst='-')
		tmp_delta=0
		if seq is not src_seq:
			for a,b in zip(src_seq,seq):
				if a != b:
					self.delta+=1
					tmp_delta+=1
					#print('N',end='')
		if tmp_delta > 0:
			if self._v: 
				print('N'*tmp_delta)
			self.changed+=1
			return tmp_rec._replace(aligned_seq_data=seq)
		return tmp_rec
	def trim_all(self, seqs):
		for s in seqs:
			yield self.trim_seq(s)
		if self._v:
			print(self.file_name,"delta",self.delta,"bp")
			print(self.file_name,"changed",self.changed,"sequences")


def write_to_file(recs, file_path,file_type="fasta",strip_aln_uid=True):
	# from itertools import imap
	if strip_aln_uid:
		seqs=map(lambda r: SeqRecord.SeqRecord(id='|'.join([r.aln_uid]+[ p for p in r.s_info if not p.startswith('ALNUID')]),
												description='',seq=r.aligned_seq_data),
				  recs)
	else:
		seqs=map(lambda r: SeqRecord.SeqRecord(id='|'.join([r.aln_uid]+r.s_info),description='',
												seq=r.aligned_seq_data),
				  recs)
	if file_type != "fasta":
		SeqIO.write(seqs,file_path,file_type)
		return
	from itertools import zip_longest
	with open(file_path,'wb') as h:
		for s in seqs:
			h.write('>%s\n' % s.id)
			c = iter(s.seq)
			for l in zip_longest(*[c]*70,fillvalue=''):
				h.write(*l,sep='')


def test_main2():
	fList=["aln-20161112PASS.old/gblocks_uce-11.fasta","aln-20161112PASS.old/gblocks_uce-11.fasta_out2"]
	outList=["trim_gblocks_uce-11.fasta","trim_gblocks_uce-11.fasta_out2"]
	worker=LocalFilter(True)
	write_to_file(worker.get_seqs(fList[0]),outList[0])
	write_to_file(worker.get_seqs(fList[1]),outList[1])

def test_main():
	fList=["aln-20161112PASS.old/gblocks_uce-11.fasta","aln-20161112PASS.old/gblocks_uce-11.fasta_out2"]
	# for testing
	from sites_local import SiteCounter, read_sequences_alt
	proc=SiteCounter()
	proc.scan_file(fList[0],"fasta")
	print("Sites Scanned")
	import sys
	sys.stdout.flush()
	delta=0
	changed=0
	for tmp_rec in read_sequences_alt(fList[0],"fasta"):
		src_seq=tmp_rec.aligned_seq_data
		foo=ambi.search(str(src_seq))
		if foo is not None:
#             print(tmp_rec.aln_uid, "Removing ambiguous")
			seq=Seq.Seq(sub_ambi(foo),src_seq.alphabet)
		else:
			seq=src_seq
		if proc.has_private_alleles(tmp_rec.aln_uid,4):
#             print(tmp_rec.aln_uid, "Scan Private")
			seq,drop=sub_priv(src_seq,
								 proc.get_private_alleles(tmp_rec.aln_uid),seq)#,
#                        tmp_rec.aligned_seq_data.alphabet)
			if len(drop) >0:
				if min(drop)>1 and get_acgt_count(src_seq[0:min(drop)-1])<10:
					seq=subst_N_sub(seq,endpos=min(drop)-1,mutable=True)
				if max(drop)<len(seq) and get_acgt_count(src_seq[max(drop):])<10:
					seq=subst_N_sub(seq,pos=max(drop),mutable=True)
		tmp_delta=0
		if seq is not src_seq:
			for a,b in zip(src_seq,seq):
				if a != b:
					delta+=1
					tmp_delta+=1
					#print('N',end='')
		if tmp_delta > 0:
			print('N'*tmp_delta)
			sys.stdout.flush()
			changed+=1
#         else:
#             print(tmp_rec.aln_uid, "Unchanged")
	print(fList[0],"delta",delta,"bp")
	print(fList[0],"changed",changed,"sequences")
	proc.scan_file(fList[1],"fasta")
	print("Sites Scanned")
	sys.stdout.flush()
	delta=0
	changed=0
	for tmp_rec in read_sequences_alt(fList[1],"fasta"):
		foo=ambi.search(str(tmp_rec.aligned_seq_data))
		if foo is not None:
#             print(tmp_rec.aln_uid, "Removing ambiguous")
			seq=Seq.Seq(sub_ambi(foo),tmp_rec.aligned_seq_data.alphabet)
		else:
			seq=tmp_rec.aligned_seq_data
		if proc.has_private_alleles(tmp_rec.aln_uid,4):
#             print(tmp_rec.aln_uid, "Scan Private")
			seq,drop=sub_priv(tmp_rec.aligned_seq_data,
								 proc.get_private_alleles(tmp_rec.aln_uid),seq)#,
#                        tmp_rec.aligned_seq_data.alphabet)
		tmp_delta=0
		if seq is not tmp_rec.aligned_seq_data:
			for a,b in zip(tmp_rec.aligned_seq_data,seq):
				if a != b:
					delta+=1
					tmp_delta+=1
					#print('N',end='')
		if tmp_delta > 0:
			print('N'*tmp_delta)
			sys.stdout.flush()
			changed+=1
#         else:
#             print(tmp_rec.aln_uid, "Unchanged")
	print(fList[1],"delta",delta,"bp")
	print(fList[1],"changed",changed,"sequences")
#     for tmp_rec in SeqIO.parse(fList[1],"fasta"):
#         foo=ambi.search(str(tmp_rec.seq))
#         if foo is not None:
#             print(tmp_rec.id)
#             print(sub_ambi(foo))


def get_counts(rec_seq):
	bar=dict(Counter(rec_seq))
	return((sum((v for k,v in bar.items() if k in 'ACGT')),
			sum((v for k,v in bar.items() if k in 'MRWSKY'))))

def get_acgt_count(rec_seq):
	return(sum((v for k,v in Counter(rec_seq).items() if k in 'ACGT')))

def subst_N_all(seq, subst='N'):
	if not isinstance(seq,str):
		seq=str(seq)
	return nongap.sub(subst, seq)

def subst_N_sub(seq, pos=None, endpos=None, mutable=False, subst='N'):
	if mutable:
		if seq.__class__=="Bio.Seq.MutableSeq":
			pass
		elif seq.__class__=="Bio.Seq.Seq":
			import Bio.Seq.MutableSeq
			seq=Bio.Seq.MutableSeq(seq,seq.alphabet)
		elif isinstance(seq,bytearray):
			pass
		else:
			mutable=False
	if pos is not None:
		if endpos is not None:
			subseq=seq[pos:endpos]
			if mutable:
				seq[pos:endpos]=nongap.sub(subst,str(subseq))
				return(seq)
		else:
			subseq=seq[pos:]
			if mutable:
				seq[pos:]=nongap.sub(subst,str(subseq))
				return(seq)
	else:
		if endpos is None:
			if mutable:
				seq[:]=nongap.sub(subst,str(subseq))
				return(seq)
			return subst_N_all(seq, subst=subst)
		subseq=seq[:endpos]
	if not isinstance(seq,str):
		subseq=str(subseq)
	res=nongap.sub(subst, subseq)
	if pos is not None:
		res = seq[:pos]+res
	if endpos is not None:
		res = res + seq[endpos:]
	return res

def sub_ambi(foo,l_sect='',r_sect='',subst='N'):
	good,poor=get_counts(foo.group(0))
	if good <= 19* poor:
		#nongap.sub('N',foo.group(0))
		if l_sect is not None:
			l_sect+=foo.string[0:foo.start()]
			if get_acgt_count(l_sect)<10:
				l_sect=nongap.sub(subst,l_sect)
		else:
			l_sect=foo.string[0:foo.start()]
		if r_sect is not None:
			r_sect=foo.string[foo.end():]+r_sect
			if get_acgt_count(r_sect)<10:
				r_sect=nongap.sub(subst,r_sect)
		else:
			r_sect=foo.string[foo.end():]
		return(''.join((l_sect,
						nongap.sub(subst,foo.group(0)),
						r_sect)))
	elif poor<4:
		return(foo.string)
	parM=foo.group(0)
	middle=acgt20.search(parM)
	if middle is None:
		return(foo.string)
	res=list()
	parL=foo.string[0:foo.start()]
	seqL=parM[0:middle.start()]
	l_srch=ambi.search(seqL)
	if l_srch is None:
		if l_sect is not None:
			res.append(l_sect)
		res.append(parL)
		res.append(seqL)
	elif l_sect is None:
		res.append(parL)
		res.append(sub_ambi(l_srch,None,None,subst=subst))
	else:
		l_sect+=parL
		if get_acgt_count(l_sect)<10:
			res.append(sub_ambi(l_srch,l_sect,None,subst=subst))
		else:
			res.append(l_sect)
			res.append(sub_ambi(l_srch,None,None,subst=subst))
	res.append(middle.group(0))
	parR=foo.string[foo.end():]
	seqR=parM[middle.end():]
	r_srch=ambi.search(seqR)
	if r_srch is None:
		res.append(seqR)
		res.append(parR)
		if r_sect is not None:
			res.append(r_sect)
	elif r_sect is None:
		res.append(sub_ambi(r_srch,None,None,subst=subst))
		res.append(parR)
	else:
		r_sect=parR+r_sect
		if get_acgt_count(r_sect)<10:
			res.append(sub_ambi(r_srch,None,r_sect,subst=subst))
		else:
			res.append(sub_ambi(r_srch,None,None,subst=subst))
			res.append(r_sect)
	return(''.join(res))

def sub_priv(seq,priv_sites,edit_seq=None, subst='-'):
	if edit_seq is None:
		edit_seq=seq
	if len(priv_sites)<4:
		return(edit_seq,list())
	priv_max=max(priv_sites)
	priv_min=min(priv_sites)-1
	if 4*len(priv_sites)>=get_acgt_count(seq[priv_min:priv_max]):
		return((edit_seq[0:priv_min]+
						nongap.sub(subst, str(edit_seq[priv_min:priv_max]))+
						edit_seq[priv_max:],priv_sites))
	#keep_sites=list()
	cur_sites=list()
	segs=list()
	cur_l =cur_r =0;
	while len(priv_sites)>3:
		while get_acgt_count(seq[(priv_sites[0]-1):priv_sites[3]])>16:
			priv_sites.pop(0)
			if len(priv_sites)<4:
				if cur_r == 0:
					return((edit_seq,list()))
				seg_iter=iter(segs)
				seg=next(seg_iter)
				for i in seg_iter:
					seg+=i
				return((seg+edit_seq[cur_r:],cur_sites))
		cur_l=priv_sites[0]-1
		segs.append(seq[cur_r:cur_l])
		for _ in range(3):
			cur_sites.append(priv_sites.pop(0))
	if len(priv_sites) > 0:
		while get_acgt_count(seq[cur_l:priv_sites[0]])<=4+4*len(cur_sites):
			cur_sites.append(priv_sites.pop(0))
			if len(priv_sites)<1:
				break
		cur_r=cur_sites[-1]
		segs.append(nongap.sub(subst, str(edit_seq[cur_l:cur_r])))
	seg_iter=iter(segs)
	seg=next(seg_iter)
	for i in seg_iter:
		seg+=i
	return(seg+edit_seq[cur_r:],cur_sites)


if __name__ == '__main__':
	main()
