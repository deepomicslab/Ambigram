import bpsmap
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd

def map_bps_junc(junc, bps_map):
    for i in junc.index:
        # print(junc.loc[i])
        junc.at[i, 'pos_5p'] = bps_map.loc[lambda df: df.chrom == junc.loc[i, 'chrom_5p']]\
                                    .loc[lambda df: df.before == junc.loc[i, 'pos_5p']]\
                                    .iloc[0].after
        junc.at[i, 'pos_3p'] = bps_map.loc[lambda df: df.chrom == junc.loc[i, 'chrom_3p']]\
                                    .loc[lambda df: df.before == junc.loc[i, 'pos_3p']]\
                                    .iloc[0].after

def map_bps_chrom_infos(chrom_infos, bps_map):
    for i in chrom_infos.index:
        chrom_infos.at[i, 'start'] = bps_map.loc[lambda df: df.chrom == chrom_infos.loc[i, 'chrom']]\
                                            .loc[lambda df: df.before == chrom_infos.loc[i, 'start']]\
                                            .iloc[0].after
        chrom_infos.at[i, 'end'] = bps_map.loc[lambda df: df.chrom == chrom_infos.loc[i, 'chrom']]\
                                          .loc[lambda df: df.before == chrom_infos.loc[i, 'end']]\
                                          .iloc[0].after

def map_bps_sv(sv, bps_map):
    for i in sv.index:
        # print(bps_map.loc[lambda df: df.before == sv.loc[i, 'pos_5p']].after)
        sv.at[i, 'pos_5p'] = bps_map.loc[lambda df: df.chrom == sv.loc[i, 'chrom_5p']]\
                                    .loc[lambda df: df.before == sv.loc[i, 'pos_5p']]\
                                    .iloc[0].after
        sv.at[i, 'pos_3p'] = bps_map.loc[lambda df: df.chrom == sv.loc[i, 'chrom_3p']]\
                                    .loc[lambda df: df.before == sv.loc[i, 'pos_3p']]\
                                    .iloc[0].after

def dedup(sv):
    sv = sv.sort_values(by=sv.columns[:6].tolist() + ['junc_reads'])
    return sv[~sv.duplicated(sv.columns[:6], keep='last')]

def segmentation(sv, chrom, start, end, id_start=1, drop_imprecise=True, drop_insertions=True):
    sv_5p = bpsmap.get_precise_sv(sv, chrom_5p=chrom, drop_imprecise=drop_imprecise, drop_insertions=drop_insertions)
    sv_3p = bpsmap.get_precise_sv(sv, chrom_3p=chrom, drop_imprecise=drop_imprecise, drop_insertions=drop_insertions)
    bps = sorted(set(bpsmap.get_breakpoints(sv_5p, sv_3p) + [start, end]))

    # bps = bpsmap.get_breakpoints(sv) + [start, end]
    # bps = list(sorted(set(bps)))
    segs = []
    for p in bps[1:-1]:
        segs.append((id_start, chrom, start, p))
        start = p
        id_start += 1
    segs.append((id_start, chrom, start, end))
    return pd.DataFrame(segs, columns=['ID', 'chrom', 'start', 'end']), id_start + 1

def update_junc_db_by_sv(sv, junc_db):
    # avg_depth = get_avg_depth(depth_tabix, chrom, start, end)
    for row in sv.itertuples():
        if row.inner_ins != '.':
            continue
        idx1 = junc_db.loc[lambda r: r.chrom_5p == row.chrom_5p]\
                      .loc[lambda r: r.pos_5p == row.pos_5p]\
                      .loc[lambda r: r.strand_5p == row.strand_5p]\
                      .loc[lambda r: r.chrom_3p == row.chrom_3p]\
                      .loc[lambda r: r.pos_3p == row.pos_3p]\
                      .loc[lambda r: r.strand_3p == row.strand_3p].index
        # print(idx)
        if idx1.empty:
            # junc_db = junc_db.append({'chrom_5p': row.chrom_5p,
            #                           'pos_5p': row.pos_5p,
            #                           'strand_5p': row.strand_5p,
            #                           'chrom_3p': row.chrom_3p,
            #                           'pos_3p': row.pos_3p,
            #                           'strand_3p': row.strand_3p,
            #                           'count': row.junc_reads / avg_depth}, ignore_index=True)
            if row.junc_reads > 5:
                junc_db = junc_db.append({'chrom_5p': row.chrom_5p,
                                          'pos_5p': row.pos_5p,
                                          'strand_5p': row.strand_5p,
                                          'chrom_3p': row.chrom_3p,
                                          'pos_3p': row.pos_3p,
                                          'strand_3p': row.strand_3p,
                                          'count': 1}, ignore_index=True)
        else:
            # junc_db.at[idx1[0], 'count'] += row.junc_reads / avg_depth
            if row.junc_reads > 5:
                junc_db.at[idx1[0], 'count'] += 1

        # idx2 = junc_db.loc[lambda r: r.chrom_5p == row.chrom_3p]\
        #               .loc[lambda r: r.pos_5p == row.pos_3p]\
        #               .loc[lambda r: r.strand_5p == '-' if row.strand_3p == '+' else '+']\
        #               .loc[lambda r: r.chrom_3p == row.chrom_5p]\
        #               .loc[lambda r: r.pos_3p == row.pos_5p]\
        #               .loc[lambda r: r.strand_3p == '-' if row.strand_5p == '+' else '+'].index
        # if idx2.empty:
        #     junc_db = junc_db.append({'chrom_5p': row.chrom_3p,
        #                               'pos_5p': row.pos_3p,
        #                               'strand_5p': '-' if row.strand_3p == '+' else '+',
        #                               'chrom_3p': row.chrom_5p,
        #                               'pos_3p': row.pos_5p,
        #                               'strand_3p': '-' if row.strand_5p == '+' else '+',
        #                               'count': 1}, ignore_index=True)
        # else:
        #     junc_db.at[idx2[0], 'count'] += 1
    return junc_db

def get_normal_junc_read_num(bam, chrom, pos, ext=5):
    n = 0
    # print(pos, ext)
    for r in bam.fetch(chrom, pos - 1, pos):
        overlapped = r.get_overlap(max(0, pos - 1 - ext), pos + ext)
        # print(r.qname, pos, ext, overlapped, pos + ext - (pos - 1 - ext))
        if overlapped == pos + ext - (pos - 1 - ext):
            n += 1
    return n

def update_junc_db_by_seg_in_chrom(segs, junc_db, bam, ext):
    # avg_depth = get_avg_depth(depth_tabix, chrom, start, end)
    for row in segs.iloc[:-1, :].itertuples():
        # print(f'Updating junc at {row.ID} {row.start} {row.end}')
        idx1 = junc_db.loc[lambda r: r.chrom_5p == row.chrom]\
                      .loc[lambda r: r.pos_5p == row.end]\
                      .loc[lambda r: r.strand_5p == '+']\
                      .loc[lambda r: r.chrom_3p == row.chrom]\
                      .loc[lambda r: r.pos_3p == row.end]\
                      .loc[lambda r: r.strand_3p == '+'].index
        # print(idx)
        if idx1.empty:
            # junc_db = junc_db.append({'chrom_5p': row.chrom,
            #                           'pos_5p': row.end,
            #                           'strand_5p': '+',
            #                           'chrom_3p': row.chrom,
            #                           'pos_3p': row.end,
            #                           'strand_3p': '+',
            #                           'count': get_normal_junc_read_num(bam, row.chrom, row.end, ext=ext) / avg_depth}, ignore_index=True)

            if get_normal_junc_read_num(bam, row.chrom, row.end, ext=ext) > 5:
                junc_db = junc_db.append({'chrom_5p': row.chrom,
                                          'pos_5p': row.end,
                                          'strand_5p': '+',
                                          'chrom_3p': row.chrom,
                                          'pos_3p': row.end,
                                          'strand_3p': '+',
                                          'count': 1}, ignore_index=True)
        else:
            # junc_db.at[idx1[0], 'count'] += get_normal_junc_read_num(bam, row.chrom, row.end, ext=ext) / avg_depth
            if get_normal_junc_read_num(bam, row.chrom, row.end, ext=ext) > 5:
                junc_db.at[idx1[0], 'count'] += 1

        # idx2 = junc_db.loc[lambda r: r.chrom_5p == row.chrom]\
        #               .loc[lambda r: r.pos_5p == row.end]\
        #               .loc[lambda r: r.strand_5p == '-']\
        #               .loc[lambda r: r.chrom_3p == row.chrom]\
        #               .loc[lambda r: r.pos_3p == row.end]\
        #               .loc[lambda r: r.strand_3p == '-'].index
        # if idx2.empty:
        #     junc_db = junc_db.append({'chrom_5p': row.chrom,
        #                               'pos_5p': row.end,
        #                               'strand_5p': '-',
        #                               'chrom_3p': row.chrom,
        #                               'pos_3p': row.end,
        #                               'strand_3p': '-',
        #                               'count': 1}, ignore_index=True)
        # else:
        #     junc_db.at[idx2[0], 'count'] += 1
    return junc_db

def write_junc_db(filename, junc_db):
    junc_db.sort_values(by=['chrom_5p', 'pos_5p', 'strand_5p', 'count']).to_csv(filename, sep='\t', index=False)

def get_avg_depth(depth, chrom, start, end):
    return sum(map(lambda x : int(x.split('\t')[-1]), depth.fetch(chrom, start, end))) / (end - start + 1)

def generate_config(filename, samplename, sv, segs, depth_tabix, bam, ext, ploidy):
    output = []
    total_depth = 0
    total_length = 0
    # for seg in segs.itertuples():
    #     total_length += seg.end - seg.start + 1
    #     total_depth += get_avg_depth(depth_tabix, seg.chrom, seg.start, seg.end) * (seg.end - seg.start + 1)
    # avg_depth = total_depth / total_length
    with open(filename, 'w') as fout:
        # fout.write(f'SAMPLE {samplename}\n')
        # fout.write(f'AVG_SEG_DP {avg_depth}\n')
        # fout.write(f'PURITY 1\n')
        # fout.write(f'AVG_PLOIDY {ploidy}\n')
        # fout.write(f'PLOIDY {ploidy}m1\n')
        # fout.write(f'SOURCE H:1\n')
        # fout.write(f'SINK H:{segs.iloc[-1].ID}\n')

        output_segs = []
        for seg in segs.itertuples():
            # print(f'Write seg {seg.ID}')
            total_length += seg.end - seg.start + 1
            seg_depth = get_avg_depth(depth_tabix, seg.chrom, seg.start, seg.end)
            total_depth += seg_depth * (seg.end - seg.start + 1)
            output_segs.append(f'SEG H:{seg.ID}:{seg.chrom}:{seg.start}:{seg.end} {seg_depth} -1')
            # fout.write(f'SEG H:{seg.ID}:{seg.chrom}:{seg.start}:{seg.end} {get_avg_depth(depth_tabix, seg.chrom, seg.start, seg.end)} -1\n')

        ins_id = len(segs) + 1
        ins_segs = []

        output_juncs = []
        juncs_depth = []
        left = next(segs.itertuples())
        for right in segs.iloc[1:].itertuples():
            support = get_normal_junc_read_num(bam, left.chrom, left.end, ext=ext)
            if support > 5:
                juncs_depth.append(support)
                output_juncs.append(f'JUNC H:{left.ID}:+ H:{right.ID}:+ {support} -1 U B')
                # fout.write(f'JUNC H:{left.ID}:+ H:{right.ID}:+ {support} -1 U B\n')
            left = right
        for row in sv.itertuples():
            # # print(row)
            # if row.junc_reads <= 5:
            #     print(row)
            #     continue
            if row.strand_5p == '+':
                if row.strand_3p == row.strand_5p:
                    left = segs.loc[lambda r: r.chrom == row.chrom_5p]\
                               .loc[lambda r: r.end == row.pos_5p]
                    right = segs.loc[lambda r: r.chrom == row.chrom_3p]\
                               .loc[lambda r: r.start == row.pos_3p]
                else:
                    left = segs.loc[lambda r: r.chrom == row.chrom_5p]\
                               .loc[lambda r: r.end == row.pos_5p]
                    right = segs.loc[lambda r: r.chrom == row.chrom_3p]\
                               .loc[lambda r: r.end == row.pos_3p]
            else:
                if row.strand_3p == row.strand_5p:
                    left = segs.loc[lambda r: r.chrom == row.chrom_5p]\
                               .loc[lambda r: r.start == row.pos_5p]
                    right = segs.loc[lambda r: r.chrom == row.chrom_3p]\
                               .loc[lambda r: r.end == row.pos_3p]
                else:
                    left = segs.loc[lambda r: r.chrom == row.chrom_5p]\
                               .loc[lambda r: r.start == row.pos_5p]
                    right = segs.loc[lambda r: r.chrom == row.chrom_3p]\
                               .loc[lambda r: r.start == row.pos_3p]
            # print(f'JUNC H:{left.ID.values[0]}:{row.strand_5p} H:{right.ID.values[0]}:{row.strand_3p} {row.junc_reads} -1 U B')
            juncs_depth.append(row.junc_reads)

            if row.inner_ins == '.':
                output_juncs.append(f'JUNC H:{left.ID.values[0]}:{row.strand_5p} H:{right.ID.values[0]}:{row.strand_3p} {row.junc_reads} -1 U B')
            else:
                ins_segs.append((ins_id, f'Ins_{ins_id}', 1, len(row.inner_ins), row.inner_ins))
                output_segs.append(f'SEG H:{ins_id}:Ins_{ins_id}:1:{len(row.inner_ins)} 1 -1')
                output_juncs.append(f'JUNC H:{left.ID.values[0]}:{row.strand_5p} H:{ins_id}:+ {row.junc_reads} -1 U B')
                output_juncs.append(f'JUNC H:{ins_id}:+ H:{right.ID.values[0]}:{row.strand_3p} {row.junc_reads} -1 U B')
                ins_id += 1
            # fout.write(f'JUNC H:{left.ID.values[0]}:{row.strand_5p} H:{right.ID.values[0]}:{row.strand_3p} {row.junc_reads} -1 U B\n')

#         print(f'SAMPLE {samplename}\n')
#         print(f'AVG_SEG_DP {total_depth * 1.0 / total_length}\n')
#         print(f'AVG_JUNC_DP {np.mean(juncs_depth)}\n')
#         print(f'PURITY 1\n')
#         print(f'AVG_PLOIDY {ploidy}\n')
#         print(f'PLOIDY {ploidy}m1\n')
#         print(f'SOURCE H:1\n')
#         print(f'SINK H:{segs.iloc[-1].ID}\n')
#         print('\n'.join(output_segs + output_juncs))

        fout.write(f'SAMPLE {samplename}\n')
        fout.write(f'AVG_SEG_DP {total_depth * 1.0 / total_length}\n')
        fout.write(f'AVG_JUNC_DP {np.mean(juncs_depth)}\n')
        fout.write(f'PURITY 1\n')
        fout.write(f'AVG_PLOIDY {ploidy}\n')
        fout.write(f'PLOIDY {ploidy}m1\n')
        fout.write(f'SOURCE H:1\n')
        fout.write(f'SINK H:{segs.iloc[-1].ID}\n')
        fout.write('\n'.join(output_segs + output_juncs) + '\n')

    if len(ins_segs) > 0:
        pd.DataFrame(ins_segs, columns=['ID', 'chrom', 'start', 'end', 'seq']).to_csv(os.path.dirname(filename) + '/' + samplename + '.inner_ins', index=False, sep='\t')
