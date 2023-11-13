#!/usr/bin/env python
from __future__ import print_function

import itertools
import re
from collections import namedtuple, Counter
from recordclass import recordclass
from Bio import SeqIO

# >>> feeder=tyranniSQL.get_connection_from_args(GD_ARGS)
# Pipeline connection made to tyranni_pipeline.
# >>> c=feeder.named_cursor()
# >>> c.execute("SELECT file_uid,len from alignment_length left join 
# alignment_sequence_ACGT_counts a using (file_uid) where a.file_uid is NULL limit 5") 



# seq allelles (14999527, 39540, 127, u'T')
SeqAlleles = namedtuple("SeqAlleles",
                        ["aln_uid", "file_uid", "site_num", "allelle"])
SeqRow = namedtuple("SeqRow",
                    ["aln_uid", "file_uid", "aligned_seq_data"])

create_tables = """CREATE TABLE private_allelles
(file_uid int, site_num int, allelle char, aln_uid int);

"""
local_columns=dict()
local_columns['file_uid']='int'
local_columns['aln_uid']='int'
local_columns['seq_uid']='int'
local_columns['site_num']='int'
local_columns['allelle']='string'
local_columns['s1']='string'
local_columns['seq_name']='string'
local_columns['snp']='string'
local_columns['site_count']='int'
local_columns['snp_count']='int'
local_columns['num_private']='int'
local_columns['ACGT_total']='int'
local_columns['A']='int'
local_columns['C']='int'
local_columns['G']='int'
local_columns['T']='int'
local_tables=dict()
local_tables['private_allelles'] = ('file_uid','site_num','allelle','aln_uid')
local_tables['alignment_site_ACGT_counts'] = ('file_uid','site_num','s1','site_count')
local_tables['alignment_depth_cache'] = ('file_uid','site_num','snp',
    'snp_count','ACGT_total','A','C','G','T','num_private')
local_tables['alignment_sequence_ACGT_counts']=('aln_uid','seq_uid','file_uid',
    'ACGT_total','A','C','G','T','num_private')
local_tables['seq_names']=('seq_uid','seq_name')



ins_priv_data = """INSERT INTO private_allelles
(file_uid,site_num,allelle,aln_uid)
VALUES (%(file_uid)s,%(site_num)s,%(allelle)s,%(aln_uid)s)"""

update_priv_count = """UPDATE alignment_sequence_ACGT_counts
SET num_private = num_private + 1
WHERE file_uid = %(file_uid)s and aln_uid = %(aln_uid)s"""

ins_allele_count = """INSERT INTO alignment_site_ACGT_counts
(file_uid,site_num,s1,site_count)
VALUES (%(file_uid)s,%(site_num)s,%(allelle)s,%(site_count)s)"""

ins_site_count = """INSERT INTO alignment_depth_cache
(file_uid,site_num,snp,snp_count,ACGT_total,A,C,G,T,num_private)
VALUES (%(file_uid)s,%(site_num)s,%(snp)s,%(snp_count)s,
%(ACGT_total)s,%(A)s,%(C)s,%(G)s,%(T)s,%(num_private)s) """

ins_seq_count = """INSERT INTO alignment_sequence_ACGT_counts
(aln_uid,seq_uid,file_uid,ACGT_total,A,C,G,T,num_private)
VALUES (%(aln_uid)s,%(seq_uid)s,%(file_uid)s,
%(ACGT_total)s,%(A)s,%(C)s,%(G)s,%(T)s,%(num_private)s)"""

select_seq_data = """SELECT aln_uid,seq_uid,file_uid,aligned_seq_data
FROM aligned_sequences WHERE file_uid = %s"""

# select_seq_data = """SELECT {0}
# FROM aligned_sequences WHERE {1} = %s
# """.format(','.join("aln_uid","seq_uid","file_uid","aligned_seq_data"),"file_uid")

select_site_data = """SELECT aln_uid,file_uid,aligned_seq_data
FROM aligned_sequences WHERE file_uid = %s"""

SeqData = namedtuple("SeqData", 
                     ["aln_uid", "seq_uid", "file_uid", 
                      "sample_name","sample_uid", "locus_uid", "aligned_seq_data"])

SeqData2 = namedtuple("SeqData2", 
                      ["aln_uid", "s_info", "file_uid", "aligned_seq_data"])

row_header=""">{SampleName}|ALNUID_{aln_uid}|SEQUID_{seq_uid}|SUID_{sample_uid}|LUID_{locus_uid}"""

max_aln_uid=9e6
max_file_uid=1e5

def get_aln_uid(aln_uid=None):
    global max_aln_uid
    if aln_uid is None:
        aln_uid=max_aln_uid
        max_aln_uid+=1
        return aln_uid
    if aln_uid >= max_aln_uid:
        max_aln_uid=aln_uid+100
    return aln_uid

def get_file_uid(file_uid=None):
    global max_file_uid
    if file_uid is None:
        file_uid=max_file_uid
        max_file_uid+=1
        return file_uid
    if file_uid >= max_file_uid:
        max_file_uid=file_uid+100
    return file_uid

def read_sequences(src_file,file_type="fasta",file_uid=None):
    seq_data=SeqIO.parse(src_file,file_type)
    #seq_cursor = sql_db.cursor()
    if file_uid is None:
        # Build a new set of records?
        file_uid=get_file_uid()
    for rec in seq_data:
        parts = rec.id.split('|')
        s_name=parts.pop(0)
        seq_id=dict(( p.split('_') for p in parts))
        yield SeqData(seq_id['ALNUID'],seq_id['SEQUID'],file_uid,
                      s_name,seq_id['SUID'],seq_id['SUID'],rec.seq)

def read_sequences_alt(src_file,file_type="fasta"):
    seq_data=SeqIO.parse(src_file,file_type)
    #seq_cursor = sql_db.cursor()
    for rec in seq_data:
        parts = rec.id.split('|')
        s_name=parts.pop(0)
        yield SeqData2(s_name, parts, 1, rec.seq)

# select_site_data = """SELECT {0}
# FROM aligned_sequences WHERE {1} = %s
# """.format(','.join("aln_uid","file_uid","aligned_seq_data"),"file_uid")

# aln_seq_ACGT_cols=("aln_uid","seq_uid","file_uid","ACGT_total","A","C","G","T","num_private")
# aln_seq_ACGT_key=("aln_uid","seq_uid","file_uid")
# aln_seq_cols=("aln_uid","seq_uid","file_uid","aligned_seq_data")
# aln_seq_ACGT_table="alignment_sequence_ACGT_counts"
# struct_insert="""INSERT INTO {0} ({1}) SELECT {1} FROM aligned_sequences WHERE {2} = %s"""

SiteAlleleData = recordclass("SiteAlleleData",
                            ["file_uid", "site_num", "A", "C", "G", "T"])
AlleleData = recordclass("AlleleData",
                        ["file_uid", "site_num", "allelle",
                         "aln_uid", "site_count"])

SiteDepthData = namedtuple("SiteDepthData", ["file_uid", "site_num", "snp",
                                             "snp_count", "ACGT_total",
                                             "A", "C", "G", "T",
                                             "num_private"])

SeqCountData = namedtuple("SeqCountData", ["aln_uid", "seq_uid", "file_uid",
                                           "ACGT_total","A","C","G","T",
                                           "num_private"])

re_not_ACGT=re.compile('[^ACGT]+')
reACGT = re.compile('[AGCT]')

def rowFilterACGT(x):
    return(reACGT.match(x.allelle) is not None)

def seq_count(row):
    seq_ACGT={'A':0,'C':0,'G':0,'T':0}
    seq_ACGT.update(Counter(re_not_ACGT.sub(u'', str(row.aligned_seq_data))))
    return(SeqCountData(row.aln_uid, row.seq_uid, row.file_uid,
                        sum(seq_ACGT.values()), num_private=0, **seq_ACGT))

def split_sites(row):
    sitesA = zip(itertools.repeat(row.aln_uid),
                            itertools.repeat(row.file_uid),
                            itertools.count(1),
                            str(row.aligned_seq_data))
    sites = itertools.starmap(SeqAlleles, sitesA)
    return(filter(rowFilterACGT, sites))

def split_sites_B(row):
    sitesA = zip(itertools.repeat(row.aln_uid),
                            itertools.repeat(row.file_uid),
                            itertools.count(1),
                            str(row.aligned_seq_data))
    sitesB = filter(lambda x: reACGT.match(x[3]) is not None, sitesA)
    sites = list(itertools.starmap(SeqAlleles, sitesB))
    return(sites)

def build_site_counts(file_uid,aln_len,rows):
    if aln_len > 0:
        res = AlnSitesV1(file_uid,aln_len)
    else:
        res = AlnSitesV2(file_uid)
    sites = itertools.chain.from_iterable(map(split_sites,rows))
    res.addmany(sites)
    return(res)

def build_site_counts_blind(rows):
    r=next(rows)
    res = AlnSitesV2(r.file_uid)
    res.addmany(split_sites(r))
    sites = itertools.chain.from_iterable(map(split_sites,rows))
    res.addmany(sites)
    return(res)

# def wedge_sites(wedge, num_wedges):
#     return(lambda x: x.site_num % num_wedges == wedge)


class SiteSeqCounter(object):
    def __init__(self, sql_db):
#        object.__init__(self,**kwargs)
        self.f_uid=None
        self.aln_sites=None
        self.db=sql_db
        
    ins_seq_count = """INSERT INTO alignment_sequence_ACGT_counts
(aln_uid,seq_uid,file_uid,ACGT_total,A,C,G,T,num_private)
VALUES (%(aln_uid)s,%(seq_uid)s,%(file_uid)s,
%(ACGT_total)s,%(A)s,%(C)s,%(G)s,%(T)s,%(num_private)s)"""

    def writer(self):
        return self.db.cursor(namedtuple=True)

    def commit(self):
        self.db.commit()

    def setup(self, src_file, file_type="fasta", file_uid=None):
        if file_uid is None:
            # Build a new set of records?
            file_uid=get_file_uid()
        self.f_uid=file_uid
        self.f_type=file_type
        self.src = src_file
#         cur=self.reader()
#         cur.execute(self.select_len_data,(file_uid,))
#         r=cur.fetchall()
#         if len(r) != 1:
#             raise Exception("File alignment length result unusable",r)
        tmp=next(read_sequences(self.src, self.f_type, self.f_uid))
        self.seq_len=len(tmp.aligned_seq_data)
        
    def scan_seq_data(self):
#         cur=self.reader()
#         cur.execute(self.select_seq_data,(self.f_uid,))
        cOut=self.writer()
        for r in read_sequences(self.src, self.f_type, self.f_uid):
            cOut.execute(self.ins_seq_count, seq_count(r)._asdict())
            yield r
        #self.commit()

    def load_site_data(self):
        seqs = self.scan_seq_data()
        self.aln_sites = build_site_counts(self.f_uid,
                                           self.seq_len,
                                           seqs)

    def push_site_data(self):
        sites = self.aln_sites.iter_sites()
        s_counts = map(lambda s: s.get_depth_data()._asdict(),sites)
        cOut=self.writer()
        cOut.executemany(ins_site_count,s_counts)
        for a_count in self.aln_sites.iter_allele_counts():
            cOut.execute(ins_allele_count,a_count._asdict())
            if is_private(a_count):
                cOut.execute(ins_priv_data,a_count._asdict())
                cOut.execute(update_priv_count,a_count._asdict())
        self.commit()
    def scan_file(self, src_file, file_type="fasta", file_uid=None):
        self.setup(src_file, file_type, file_uid)
        self.load_site_data()
        self.push_site_data()
        return(self.f_uid)

class SiteCounter(object):
    def __init__(self):
#        object.__init__(self,**kwargs)
#         self.f_uid=None
        self.aln_sites=None
#         self.db=sql_db
#         
#     ins_seq_count = """INSERT INTO alignment_sequence_ACGT_counts
# (aln_uid,seq_uid,file_uid,ACGT_total,A,C,G,T,num_private)
# VALUES (%(aln_uid)s,%(seq_uid)s,%(file_uid)s,
# %(ACGT_total)s,%(A)s,%(C)s,%(G)s,%(T)s,%(num_private)s)"""
# 
#     def writer(self):
#         return self.db.cursor(namedtuple=True)
# 
#     def commit(self):
#         self.db.commit()

    def setup(self, src_file, file_type="fasta"):
        self.priv_alleles=dict()
        self.f_type=file_type
        self.src = src_file
#         cur=self.reader()
#         cur.execute(self.select_len_data,(file_uid,))
#         r=cur.fetchall()
#         if len(r) != 1:
#             raise Exception("File alignment length result unusable",r)
        tmp=next(read_sequences_alt(self.src, self.f_type))
        self.seq_len=len(tmp.aligned_seq_data)
        
    def scan_seq_data(self):
#         cur=self.reader()
#         cur.execute(self.select_seq_data,(self.f_uid,))
        for r in read_sequences_alt(self.src, self.f_type):
            self.priv_alleles[r.aln_uid]=list()
            yield r
        #self.commit()

    def load_site_data(self):
        seqs = self.scan_seq_data()
        self.aln_sites = build_site_counts(1,
                                           self.seq_len,
                                           seqs)

#     def get_site_data(self):
#         sites = self.aln_sites.iter_sites()
#         s_counts = itertools.imap(lambda s: s.get_depth_data()._asdict(),sites)
#         cOut=self.writer()
#         cOut.executemany(ins_site_count,s_counts)
#         for a_count in self.aln_sites.iter_allele_counts():
#             cOut.execute(ins_allele_count,a_count._asdict())
#             if is_private(a_count):
#                 cOut.execute(ins_priv_data,a_count._asdict())
#                 cOut.execute(update_priv_count,a_count._asdict())
#         self.commit()
    def scan_file(self, src_file, file_type="fasta"):
        self.setup(src_file, file_type)
        self.load_site_data()
        self.scan_private_alleles()

    def scan_private_alleles(self):
        for a_count in self.aln_sites.iter_allele_counts():
            if is_private(a_count):
                self.priv_alleles[a_count.aln_uid].append(a_count.site_num)
    def has_private_alleles(self,aln_uid,min_count=1):
        return( len(self.priv_alleles[aln_uid]) )
    def get_private_alleles(self,aln_uid):
        return( sorted(self.priv_alleles[aln_uid]) )




class AlleleCount(AlleleData):
    def __new__(_cls,file_uid, site_num, allelle,
                         aln_uid=0, site_count=0):
        return AlleleData.__new__(_cls,file_uid, site_num, allelle,
                         aln_uid, site_count)

    def add(self, row):
        if self.site_count == 0:
            self.aln_uid = row.aln_uid
        else:
            self.aln_uid = 0
        self.site_count += 1

    def merge(self, row):
        if self.site_count == 0:
            self.aln_uid = row.aln_uid
        else:
            self.aln_uid = 0
        self.site_count += row.site_count


def is_private(site_row):
    return site_row.site_count == 1


def not_empty(site):
    return site.site_count > 0


class SiteCounts(SiteAlleleData):
    __slots__ = ()
    def __new__(_cls,file_uid, site_num):
        return SiteAlleleData.__new__(_cls,file_uid, site_num,
                                      AlleleCount(file_uid, site_num, 'A'),
                                      AlleleCount(file_uid, site_num, 'C'),
                                      AlleleCount(file_uid, site_num, 'G'),
                                      AlleleCount(file_uid, site_num, 'T'))
#     def __init__(self, f_uid, site):
#         SiteAlleleData.__init__(self, f_uid, site,
#                                 AlleleCount(f_uid, site_num, 'A'),
#                                 AlleleCount(f_uid, site_num, 'C'),
#                                 AlleleCount(f_uid, site_num, 'G'),
#                                 AlleleCount(f_uid, site_num, 'T'))

    def add(self, row):
        if row.allelle == 'A':
            self.A.add(row)
            return
        if row.allelle == 'C':
            self.C.add(row)
            return
        if row.allelle == 'G':
            self.G.add(row)
            return
        if row.allelle == 'T':
            self.T.add(row)
            return

    def is_empty(self):
        if self.A.site_count > 0:
            return False
        if self.C.site_count > 0:
            return False
        if self.G.site_count > 0:
            return False
        if self.T.site_count > 0:
            return False
        return True

    def iter_alleles(self):
        if self.A.site_count > 0:
            yield self.A
        if self.C.site_count > 0:
            yield self.C
        if self.G.site_count > 0:
            yield self.G
        if self.T.site_count > 0:
            yield self.T

    def iter_counts(self):
        yield self.A.site_count
        yield self.C.site_count
        yield self.G.site_count
        yield self.T.site_count

    def get_depth_data(self):
        alleles = sorted(self.iter_alleles(),
                         key=lambda x: x.site_count,
                         reverse=True)
        snp = ''.join(map(lambda x: x.allelle, alleles))
        return SiteDepthData(self.file_uid, self.site_num,
                             snp, len(alleles), sum(self.iter_counts()),
                             self.A.site_count,
                             self.C.site_count,
                             self.G.site_count,
                             self.T.site_count,
                             len(filter(is_private, alleles)))


def merge_count(c1, c2):
    if c1.is_empty():
        return c2
    if c2.is_empty():
        return c1
    c1.A.merge(c2.A)
    c1.C.merge(c2.C)
    c1.G.merge(c2.G)
    c1.T.merge(c2.T)
    return c1



class AlnSitesV1(object):
    def __init__(self, f_uid, seq_len):
        self.sites = [SiteCounts(f_uid, site_i + 1)
                      for site_i in range(seq_len)]

    def iter_allele_counts(self):
        return itertools.chain.from_iterable(map(lambda x: x.iter_alleles(),
                                                            self.sites))

    def iter_sites(self):
        return iter(self.sites)

    def add(self, site_row):
        self.sites[site_row.site_num - 1].add(site_row)

    def addmany(self, site_rows):
        for r in site_rows:
            self.sites[r.site_num - 1].add(r)

    def merge(self, aln_sites):
        if not isinstance(aln_sites, AlnSitesV1):
            for s in aln_sites.iter_sites():
                self.sites[s.site_num
                           - 1] = merge_count(self.sites[s.site_num - 1], s)
            return
        junk = zip(self.sites, aln_sites.sites)
        self.sites = [merge_count(a, b) for a, b in junk]


class AlnSitesV2(object):
    def __init__(self, f_uid):
        self.sites = dict()

    def iter_allele_counts(self):
        return itertools.chain.from_iterable(self.sites.itervalues())

    def iter_sites(self):
        return self.sites.itervalues()

    def add(self, site_row):
        if site_row.site_num not in self.sites.keys():
            self.sites[site_row.site_num
                       ] = SiteCounts(site_row.file_uid,
                                      site_row.site_num)
        self.sites[site_row.site_num].add(site_row)

    def addmany(self, site_rows):
        for r in site_rows:
            if r.site_num not in self.sites.keys():
                self.sites[r.site_num
                           ] = SiteCounts(r.file_uid,
                                          r.site_num)
            self.sites[r.site_num].add(r)

    def merge(self, aln_sites):
        for s in aln_sites.iter_sites():
            if s.site_num not in self.sites.keys():
                self.sites[s.site_num] = s
            else:
                self.sites[s.site_num
                           ] = merge_count(self.sites[s.site_num],
                                           s)
            return
        junk = zip(self.sites, aln_sites.sites)
        self.sites = [merge_count(a, b) for a, b in junk]

        
def insert_aln_sites(cur,aln_sites_obj):
    pass
