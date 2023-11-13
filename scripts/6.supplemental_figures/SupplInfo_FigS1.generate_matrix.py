import re
import glob


def get_seq(aln, inds):
    f = open(aln, 'r')
    ids = {}
    for l in f:
        if re.search('>', l.strip()):
            id = re.search('>(\S+)', l).group(1)
            id = re.sub('_R_', '', id)
            if id not in inds:
                inds[id] = 1
            ids[id] = 1
    f.close()
    return inds, ids
    
dir = '/nfs/turbo/lsa-rabosky/Lab/squamate_tree/all_trim_300_0.05_0.3/'
alns = glob.glob(dir + '*aln')

inds = {}
map = {}
for aln in alns:
    loc = re.sub(dir, '', aln)
    loc = re.sub('.fasta.*', '', loc)
    inds, ids = get_seq(aln, inds)
    map[loc] = ids

out = '/nfs/turbo/lsa-rabosky/Lab/squamate_tree/sample_matrix.csv'
o = open(out, 'w')

loci = sorted(list(map.keys()))
o.write('ind,' + ','.join(loci) + '\n') 
for ind in inds:
    vals = [ind]
    for loc in loci:
        if ind in map[loc]:
            vals.append('1')
        else:
            vals.append('0')
    o.write(','.join(vals) + '\n')
o.close()
    
