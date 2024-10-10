
def refine_chr(chr_):
    chr_ = chr_.strip().lower()
    for _ in ['chromosome', 'chrom', 'chr']:
        chr_ = chr_.replace(_, '')
    return {'23': 'x', '24': 'y', '25': 'mt', 'm': 'mt'}.get(chr_) or chr_

# 8:	associated_gene-50kb	PALMD

# 9:	associated_gene-4DGenome	PALMD
# 10:	Cell_Tissue	IMR90
# 11:	Detection_method	Hi-C
# 12:	PMID	24141950

# 13:	closest_gene	PALMD
# 14:	distance	20452

def process(fn):
    res = {}  # k = region chr:s:e, v = [clist, fantom, asso]
    site_anno = {}  # key = center point, v = [fantom_v, asso_v]
    with open(fn) as f:
        f.readline()
        for i in f:
            line = i[:-1].split('\t')
            chr_, s, e, clist, fantom, asso = line[1:7]
            fd_gn, fd_cell, fd_method, fd_pmid = line[8:12]
            fd_gn = ','.join(sorted(fd_gn.split(',')))
            fd_cell = ','.join(sorted(fd_cell.split(',')))
            fd_method = ','.join(sorted(fd_method.split(',')))
            fd_pmid = ','.join(sorted(fd_pmid.split(',')))
            close_50k = ','.join(sorted(line[7].split(',')))
            closest = ','.join(sorted(line[12].split(',')))
            distance = line[13]
            fd = [fd_gn, fd_cell, fd_method, fd_pmid]
            chr_ = refine_chr(chr_)
            clist = clist.split(',')
            fantom = fantom.split(',')
            asso = asso.split(';')
            clist_with_idx = sorted(list(enumerate(clist)), key=lambda x: int(x[1]))
            clist_new, fantom_new, asso_new = [], [], []
            for idx, center in clist_with_idx:
                clist_new.append(center)
                site_anno.setdefault(chr_, {})[int(center)] = [fantom[idx], asso[idx]]
                fantom_new.append(fantom[idx])
                asso_new.append(','.join(sorted(asso[idx].split(',')))) 
            res[f'{chr_}:{s}:{e}'] = [clist_new, fantom_new, asso_new, fd, close_50k, closest, distance]
    return res, site_anno

# under bioinfo2,  /data/nrsa/testv2/encode_gtf/eRNA
fn1 = '/data/nrsa/testv1/encode/eRNA/Enhancer.txt'
fn2 = '/data/nrsa/testv2/encode_gtf/eRNA/Enhancer.txt'

res1, site1 = process(fn1) # res1 = 1084
res2, site2 = process(fn2)  # res2 = 1085
del res2['y:19441936:19445693'] # version1 missing the last item

# compare center list
diff_center = []
diff_fantom = []
diff_asso = []
diff_4d = []
diff_50k = []
diff_closest = []
diff_distance = []
for k in res1:
    v1 = res1[k]
    v2 = res2[k]
    for idx, diff_l in zip(range(7), [diff_center, diff_fantom, diff_asso, diff_4d, diff_50k, diff_closest, diff_distance]):
        if v1[idx] != v2[idx]:
            diff_l.append([k, v1[0], v1[idx], v2[idx]])

len(diff_center) # 0 , exact match
len(diff_fantom) # 0 , exact match
len(diff_asso) # 0 , exact match
len(diff_4d) # 0 , exact match
len(diff_50k) # 0 , exact match
len(diff_closest) # 0 , exact match
len(diff_distance) # 0 , exact match

import pickle
fn_fantom_sites = '/data/nrsa/fantom_sites.pkl'
fn_asso_sites = '/data/nrsa/enh_sites.pkl'

with open(fn_fantom_sites, 'rb') as f:
    fantom_sites = pickle.load(f)
with open(fn_asso_sites, 'rb') as f:
    asso_sites = pickle.load(f)

fn_4d = '/data/nrsa/code/annotation/4DGenome-hg19.txt'
d_4d = {}
if fn_4d is not None:
    # chr1	557489	560146	chr5	134284878	134293544	PCBD2;	MCF7	ChIA-PET	19890323
    with open(fn_4d) as f:
        # no header
        for i in f:
            line = i[:-1].split('\t')
            chr_ = refine_chr(line[0])
            d_4d.setdefault(chr_, []).append([line[_] for _ in [1, 2, 6, 7, 8, 9]]) #  s, e, gn, cell, method, pmid 
