import pandas as pd
import numpy as np
import pysam


from sklearn.neighbors import KDTree
from scipy.spatial.distance import pdist, squareform

def merge_sv_tgs2sgs(sgs, tgs, thres):
    res = pd.DataFrame(columns=['chrom_5p', 'pos_5p', 'strand_5p',
                                 'chrom_3p', 'pos_3p', 'strand_3p',
                                 'inner_ins', 'span_reads', 'junc_reads',
                                 'id', 'qual', 'filter', 'meta_info', 'anno_info'])
    res = res.astype({'chrom_5p': str, 'pos_5p': np.int64, 'strand_5p': str,
                 'chrom_3p': str, 'pos_3p': np.int64, 'strand_3p': str,
                 'inner_ins': str, 'span_reads': np.int64, 'junc_reads': np.int64,
                 'id': str, 'qual': np.float64, 'filter': str, 'meta_info': str, 'anno_info': str})
    n_found = 0
    n_found_comp = 0
    for row in tgs.itertuples():
#         print('----------------')
#         print(res)
#         print('##########################')
#         print(row)
        found = res.loc[lambda r: r.chrom_5p == row.chrom_5p]\
                    .loc[lambda r: r.chrom_3p == row.chrom_3p]\
                    .loc[lambda r: r.strand_5p == row.strand_5p]\
                    .loc[lambda r: r.strand_3p == row.strand_3p]
        found_comp = res.loc[lambda r: r.chrom_5p == row.chrom_3p]\
                        .loc[lambda r: r.chrom_3p == row.chrom_5p]\
                        .loc[lambda r: r.strand_5p == ('+' if row.strand_3p == '-' else '-')]\
                        .loc[lambda r: r.strand_3p == ('+' if row.strand_5p == '-' else '-')]
#         print('+++++++++++++++++++++++++')
#         print(found)
        if found.empty and found_comp.empty:
            found = sgs.loc[lambda r: r.chrom_5p == row.chrom_5p]\
                        .loc[lambda r: r.chrom_3p == row.chrom_3p]\
                        .loc[lambda r: r.strand_5p == row.strand_5p]\
                        .loc[lambda r: r.strand_3p == row.strand_3p]
            found_comp = sgs.loc[lambda r: r.chrom_5p == row.chrom_3p]\
                            .loc[lambda r: r.chrom_3p == row.chrom_5p]\
                            .loc[lambda r: r.strand_5p == ('+' if row.strand_3p == '-' else '-')]\
                            .loc[lambda r: r.strand_3p == ('+' if row.strand_5p == '-' else '-')]
#         print('|||||||||||||||||||||||')
#         print(found)
#         if not found.empty or not found_comp.empty:
#             return row, found, found_comp
        if len(found) > 0:
            n_found += 1
            dist = (found.pos_5p - row.pos_5p).abs() + (found.pos_3p - row.pos_3p).abs()
            if dist.min() <= thres * 2:
#                 print('*************')
#                 print(found)
                rf = found.loc[dist.idxmin()]
                idx1 = res.loc[lambda r: r.chrom_5p == rf.chrom_5p]\
                                 .loc[lambda r: r.pos_5p == rf.pos_5p]\
                                 .loc[lambda r: r.strand_5p == rf.strand_5p]\
                                 .loc[lambda r: r.chrom_3p == rf.chrom_3p]\
                                 .loc[lambda r: r.pos_3p == rf.pos_3p]\
                                 .loc[lambda r: r.strand_3p == rf.strand_3p].index
                idx1_comp = res.loc[lambda r: r.chrom_5p == rf.chrom_3p]\
                                 .loc[lambda r: r.pos_5p == rf.pos_3p]\
                                 .loc[lambda r: r.strand_5p == ('+' if rf.strand_3p == '-' else '-')]\
                                 .loc[lambda r: r.chrom_3p == rf.chrom_5p]\
                                 .loc[lambda r: r.pos_3p == rf.pos_5p]\
                                 .loc[lambda r: r.strand_3p == ('+' if rf.strand_5p == '-' else '-')].index
                if idx1.empty and idx1_comp.empty:
                    res = res.append(rf)
                else:
                    if not idx1.empty:
                        res.at[idx1[0], 'junc_reads'] += rf.junc_reads
                    elif not idx1_comp.empty:
                        res.at[idx1_comp[0], 'junc_reads'] += rf.junc_reads
#                 print(f'sgs: {r.chrom_5p} {r.pos_5p} {r.strand_5p} {r.chrom_3p} {r.pos_3p} {r.strand_3p}')
        elif len(found_comp) > 0:
            n_found_comp += 1
            dist = (found_comp.pos_5p - row.pos_3p).abs() + (found_comp.pos_3p - row.pos_5p).abs()
            if dist.min() <= thres * 2:
                rf = found_comp.loc[dist.idxmin()]
                idx1 = res.loc[lambda r: r.chrom_5p == rf.chrom_5p]\
                                 .loc[lambda r: r.pos_5p == rf.pos_5p]\
                                 .loc[lambda r: r.strand_5p == rf.strand_5p]\
                                 .loc[lambda r: r.chrom_3p == rf.chrom_3p]\
                                 .loc[lambda r: r.pos_3p == rf.pos_3p]\
                                 .loc[lambda r: r.strand_3p == rf.strand_3p].index
                idx1_comp = res.loc[lambda r: r.chrom_5p == rf.chrom_3p]\
                                 .loc[lambda r: r.pos_5p == rf.pos_3p]\
                                 .loc[lambda r: r.strand_5p == ('+' if rf.strand_3p == '-' else '-')]\
                                 .loc[lambda r: r.chrom_3p == rf.chrom_5p]\
                                 .loc[lambda r: r.pos_3p == rf.pos_5p]\
                                 .loc[lambda r: r.strand_3p == ('+' if rf.strand_5p == '-' else '-')].index
                if idx1.empty and idx1_comp.empty:
                    res = res.append(rf)
                else:
                    if not idx1.empty:
                        res.at[idx1[0], 'junc_reads'] += rf.junc_reads
                    elif not idx1_comp.empty:
                        res.at[idx1_comp[0], 'junc_reads'] += rf.junc_reads
        else:
            rf = row
            idx1 = res.loc[lambda r: r.chrom_5p == rf.chrom_5p]\
                                 .loc[lambda r: r.pos_5p == rf.pos_5p]\
                                 .loc[lambda r: r.strand_5p == rf.strand_5p]\
                                 .loc[lambda r: r.chrom_3p == rf.chrom_3p]\
                                 .loc[lambda r: r.pos_3p == rf.pos_3p]\
                                 .loc[lambda r: r.strand_3p == rf.strand_3p].index
            idx1_comp = res.loc[lambda r: r.chrom_5p == rf.chrom_3p]\
                                 .loc[lambda r: r.pos_5p == rf.pos_3p]\
                                 .loc[lambda r: r.strand_5p == ('+' if rf.strand_3p == '-' else '-')]\
                                 .loc[lambda r: r.chrom_3p == rf.chrom_5p]\
                                 .loc[lambda r: r.pos_3p == rf.pos_5p]\
                                 .loc[lambda r: r.strand_3p == ('+' if rf.strand_5p == '-' else '-')].index
            if idx1.empty and idx1_comp.empty:
                od = rf._asdict()
                name = od['Index']
                od.pop('Index')
                ods = pd.Series(od)
                ods.name = name
                res = res.append(ods)
            else:
                if not idx1.empty:
                    res.at[idx1[0], 'junc_reads'] += rf.junc_reads
                elif not idx1_comp.empty:
                    res.at[idx1_comp[0], 'junc_reads'] += rf.junc_reads

#                 print(f'sgs_comp: {r.chrom_5p} {r.pos_5p} {r.strand_5p} {r.chrom_3p} {r.pos_3p} {r.strand_3p}')

    return n_found, n_found_comp, res

def concat_sv(sv_list_filename):
    df = pd.DataFrame()
    with open(sv_list_filename, 'r') as fin:
        for line in fin:
            df = df.append(read_sv(line[:-1]))
    return df

def read_sv(file_name):
    return pd.read_csv(file_name, header=None, sep='\t',
                          names=['chrom_5p', 'pos_5p', 'strand_5p',
                                 'chrom_3p', 'pos_3p', 'strand_3p',
                                 'inner_ins', 'span_reads', 'junc_reads',
                                 'id', 'qual', 'filter', 'meta_info', 'anno_info'])


def get_precise_sv(sv_df, chrom_5p=None, start_5p=None, end_5p=None,
                   chrom_3p=None, start_3p=None, end_3p=None,
                   drop_imprecise=True, drop_insertions=True,
                   support_thres=5):
    # depth_tabix = pysam.TabixFile(depth_filename)
    # avg_depth = get_avg_depth(depth_tabix, chrom, start, end)
    # res_df = pd.read_table(sv_filename, header=None,
    #                       names=['chrom_5p', 'pos_5p', 'strand_5p',
    #                              'chrom_3p', 'pos_3p', 'strand_3p',
    #                              'inner_ins', 'span_reads', 'junc_reads',
    #                              'id', 'qual', 'filter', 'meta_info', 'anno_info'])
    # print(next(sv_df.itertuples()).chrom_5p, chrom)
    res_df = sv_df
    if chrom_5p:
        res_df = res_df.loc[lambda row: row.chrom_5p == chrom_5p]
        if start_5p:
            res_df = res_df.loc[lambda row: row.pos_5p >= start_5p]
        if end_5p:
            res_df = res_df.loc[lambda row: row.pos_5p <= end_5p]

    if chrom_3p:
        res_df = res_df.loc[lambda row: row.chrom_3p == chrom_3p]
        if start_3p:
            res_df = res_df.loc[lambda row: row.pos_3p >= start_3p]
        if end_3p:
            res_df = res_df.loc[lambda row: row.pos_3p <= end_3p]

    res_df = res_df.loc[lambda row: row.junc_reads >= support_thres]
    if drop_imprecise:
        res_df = res_df.loc[lambda row: ~row.meta_info.str.contains('IMPRECISE')]
    if drop_insertions:
        res_df = res_df.loc[lambda row: ~row.meta_info.str.contains('INSERTION')]
    return res_df

# def get_precise_sv_seeksv(sv_filename, chrom='chr6', start=28460000, end=33500000, support_thres=5):
#     sv_df = pd.read_table(sv_filename, skiprows=1, header=None,
#                           usecols=[0, 1, 2, 3, 4, 5, 6, 7],
#                           names=['chrom_5p', 'pos_5p', 'strand_5p', 'left_read',
#                                  'chrom_3p', 'pos_3p', 'strand_3p', 'right_read'])
#     res_df = sv_df.loc[lambda row: row.chrom_5p == chrom]\
#                   .loc[lambda row: row.pos_5p >= start]\
#                   .loc[lambda row: row.pos_5p <= end]\
#                   .loc[lambda row: row.chrom_3p == chrom]\
#                   .loc[lambda row: row.pos_3p >= start]\
#                   .loc[lambda row: row.pos_3p <= end]\
#                   .loc[lambda row: row.left_read >= support_thres]\
#                   .loc[lambda row: row.right_read >= support_thres]
#     return res_df

def get_breakpoints(sv_5p, sv_3p):
    return sorted(set(sv_5p.pos_5p).union(sv_3p.pos_3p))

def get_breakpoints_from_list(sv_list_filename, chrom, start, end, support_thres=5):
    bps_set = set()
    n = 1
    with open(sv_list_filename, 'r') as fin:
        for line in fin:
            sv_filename, depth_filename = line[:-1].split()
            print(f'{n} {line}')
            n += 1
            sv = get_precise_sv(sv_filename, depth_filename, chrom, start, end, support_thres)
            bps_set = bps_set.union(get_breakpoints(sv))
    return np.array(sorted(bps_set))

def count_neighbor(arr, r=20):
	dist = squareform(pdist(arr))
	return np.array([sum(d < r) for d in dist])

def map_bps(bps, r):
    bps_rs = bps.reshape(-1, 1)
    kdt = KDTree(bps_rs)
    ns = kdt.query_radius(bps_rs, r=10)

    inters = []
    inter = set(ns[0])
    for i in range(1, len(ns)):
    	if inter.intersection(set(ns[i])):
    		inter = inter.union(set(ns[i]))
    	else:
    		inters.append(list(inter))
    		inter = set(ns[i])
    inters.append(list(inter))
    bps_map = []
    for a in inters:
    	p = bps_rs[a]
    	c = count_neighbor(p, r)
    	pivot = p[c.argmax()][0]
    	for n in p:
    		bps_map.append((n[0], pivot))
    return bps_map

def write_bps_map(filename, bps_map):
    # bps_map = [(chrom, *t) for t in bps_map]
    pd.DataFrame(bps_map, columns=['chrom', 'before', 'after'])\
      .sort_values(by=['chrom', 'before'])\
      .to_csv(filename, sep='\t', index=False)
