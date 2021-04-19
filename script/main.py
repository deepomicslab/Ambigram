import argparse
import sys, os

os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

class MainArgParser:
    def __init__(self):
        parser = argparse.ArgumentParser(prog='prelocalhap')
        parser.add_argument(dest='subfunc', help='Subcommands: ')
        args = parser.parse_args(sys.argv[1:2])
        getattr(self, args.subfunc)()

        
    def unmap2ins(self):
        import pandas as pd
        parser = argparse.ArgumentParser(description='Convert unmapped segments to insertions')
        parser.add_argument('-f', '--sv-file',
                            dest='sv_file',
                            required=True,
                            help='SV file')
        parser.add_argument('-u', '--unmapped-common-string',
                            dest='unmapped_common_str',
                            required=True,
                            help='Unmapped common string')
        parser.add_argument('-r', '--ref-common=string',
                            dest='ref_common_str',
                            required=True,
                            help='Reference common string')
        parser.add_argument('-o', '--out-sv',
                            dest='out_sv',
                            required=True,
                            help='Output path for SV.')
        args = parser.parse_args(sys.argv[2:])
        
        sv = pd.read_csv(args.sv_file, sep='\t',
           names=['chrom_5p', 'pos_5p', 'strand_5p',
                  'chrom_3p', 'pos_3p', 'strand_3p',
                  'inner_ins', 'span_reads', 'junc_reads',
                  'id', 'qual', 'filter', 'meta_info', 'anno_info'])
    
        segs = set(sv.query(f'chrom_3p.str.contains("{args.unmapped_common_str}")').chrom_3p) \
                .union(sv.query(f'chrom_5p.str.contains("{args.unmapped_common_str}")').chrom_5p)

        d = []
        for s in segs:
    #         print(f'{sample} {s}')
            df = sv.query(f'chrom_5p == "{s}" | chrom_3p == "{s}"')
            if df.chrom_5p.unique().size == 1 or df.chrom_3p.unique().size == 1:
                new = df.iloc[0]
                if df.chrom_5p.unique().size == 1:
                    new.chrom_5p = df.iloc[1].chrom_3p
                    new.pos_5p = df.iloc[1].pos_3p
                    new.strand_5p = ('+' if df.iloc[1].strand_3p == '-' else '-')
                elif df.chrom_3p.unique().size == 1:
                    new.chrom_3p = df.iloc[1].chrom_5p
                    new.pos_3p = df.iloc[1].pos_5p
                    new.strand_3p = ('+' if df.iloc[1].strand_5p == '-' else '-')
            else:
                new = df.query(f'chrom_3p == "{s}"').iloc[0]
                new.chrom_3p = df.query(f'chrom_5p == "{s}"').iloc[0].chrom_3p
                new.pos_3p = df.query(f'chrom_5p == "{s}"').iloc[0].pos_3p
                new.strand_3p = df.query(f'chrom_5p == "{s}"').iloc[0].strand_3p
            new.inner_ins = s
            new.junc_reads = df.junc_reads.min()
            d.append(new)

        d = pd.DataFrame(d)
        d = d.append(sv.query(f'chrom_3p.str.contains("{args.ref_common_str}") & chrom_5p.str.contains("{args.ref_common_str}")'))
        d.to_csv(args.out_sv, index=False, header=False, sep='\t')
    
    def bpsmap(self):
        import bpsmap
        import numpy as np
        import pandas as pd
        parser = argparse.ArgumentParser(description='Group and merge breakpoints, output mapping file')
        parser.add_argument('-l', '--sv-list',
                            dest='sv_list',
                            required=True,
                            help='SV file list')
        parser.add_argument('-r', '--region',
                            dest='region',
                            default=None,
                            help='Region')
        parser.add_argument('-S', '--chrom-sv',
                            dest='chrom_sv',
                            required=False,
                            default=None,
                            help='Chromosome SV')
        parser.add_argument('-C', '--chrom-info',
                            dest='chrom_info',
                            required=True,
                            help='Chromosome information')
        parser.add_argument('-i', '--keep-imprecise',
                            dest='keep_imprecise',
                            action='store_true',
                            default=True,
                            help='Keep imprecise SV')
        parser.add_argument('-I', '--keep-insertions',
                            dest='keep_insertions',
                            action='store_true',
                            default=True,
                            help='Keep insertions')
        parser.add_argument('-o', '--out-bps-map',
                            dest='bps_map_out',
                            required=True,
                            help='Output path for breakpoint map.')
        args = parser.parse_args(sys.argv[2:])

        if args.chrom_sv is not None:
            chrom_junc = pd.read_csv(args.chrom_sv,
                                     sep='\t',
                                     names=[
                                         'chrom_5p', 'pos_5p', 'strand_5p',
                                         'chrom_3p', 'pos_3p', 'strand_3p', 'count'
                                         ]
                                    )
        else:
            chrom_junc = None
        seek_sv = False
        if seek_sv:
            sv_df = pd.read_table(args.sv_list, skiprows=1, header=None,
                                  usecols=[0, 1, 2, 3, 4, 5, 6, 7],
                                  names=['chrom_5p', 'pos_5p', 'strand_5p', 'left_read',
                                         'chrom_3p', 'pos_3p', 'strand_3p', 'right_read'])
            # sv_df = pd.read_csv(args.sv_list)
        else:
            sv_df = bpsmap.concat_sv(args.sv_list)
        chrom_infos = pd.read_csv(args.chrom_info,dtype={'end':np.int64}, sep='\t')
        # print(chrom_junc.head())
        # print(sv_df.head())
        # print(chrom_infos.head())
        if args.region:
            ## TODO: if region not None
            pass
        else:
            chroms = list(sv_df.chrom_5p.unique()) + list(sv_df.chrom_3p.unique()) + chrom_infos.chrom.tolist()
            if chrom_junc is not None:
                chroms = sorted(set(
                    list(chrom_junc.chrom_5p.unique()) +
                    list(chrom_junc.chrom_3p.unique()) +
                    chroms))
            else:
                # print(chroms)
                print(list(sv_df.columns))

            chroms = sorted(set(chroms))
            bps_map = []
            for chrom in chroms:
                sv_5p = bpsmap.get_precise_sv(sv_df, chrom_5p=chrom, drop_imprecise=not args.keep_imprecise, drop_insertions=not args.keep_insertions)
                sv_3p = bpsmap.get_precise_sv(sv_df, chrom_3p=chrom, drop_imprecise=not args.keep_imprecise, drop_insertions=not args.keep_insertions)
                if chrom_junc is not None:
                    chrom_junc_5p = chrom_junc.loc[lambda df: df.chrom_5p == chrom]
                    chrom_junc_3p = chrom_junc.loc[lambda df: df.chrom_3p == chrom]
                print(chrom)
                chrom_info = chrom_infos.loc[lambda df: df.chrom == chrom].iloc[0]
                bps = bpsmap.get_breakpoints(sv_5p, sv_3p) + [chrom_info.start, chrom_info.end]
                if chrom_junc is not None:
                    bps = np.array(sorted(set(
                        bps +
                        bpsmap.get_breakpoints(chrom_junc_5p, chrom_junc_3p)
                        )))
                else:
                    bps = np.array(sorted(set(bps)))
                bps_map.extend([(chrom, *t) for t in bpsmap.map_bps(bps, 10)])
            bpsmap.write_bps_map(args.bps_map_out, bps_map)


        # chrom = args.region.split(':')[0]
        # start, end = [int(i) for i in args.region.split(':')[1].split('-')]
        # bps = bpsmap.get_breakpoints_from_list(args.sv_list, chrom, start, end, svaba=not args.is_seeksv)
        # bps_map = bpsmap.map_bps(bps)
        # bpsmap.write_bps_map(args.bps_map_out, bps_map, chrom)

    def config(self):
        import bpsmap
        import config
        import pysam
        import pandas as pd
        import numpy as np
        parser = argparse.ArgumentParser(description='Generate localhap config for each individual')
        parser.add_argument('-f', '--sv-file',
                            dest='sv_file',
                            required=True,
                            help='Individual SV file')
        parser.add_argument('-b', '--bam-file',
                            dest='bam_file',
                            required=True,
                            help='Individual BAM file')
        parser.add_argument('-S', '--seeksv',
                            dest='is_seeksv',
                            required=False,
                            default=False,
                            action='store_true',
                            help='Whether seeksv results')
        parser.add_argument('-m', '--bps-map',
                            dest='bps_map',
                            required=True,
                            help='Breakpoint map file')
        # group.add_argument('-J', '--create-junc-db',
        #                     dest='new_junc_db',
        #                     help='New junction database')
        parser.add_argument('-j', '--junc-db',
                            dest='junc_db',
                            required=True,
                            help='Junction database')
        parser.add_argument('-d', '--depth-tabix',
                            dest='depth_file',
                            required=True,
                            help='Tabixed depth for counting supports')
        parser.add_argument('-C', '--chrom-info',
                            dest='chrom_info',
                            required=True,
                            help='Chromosome information')
        parser.add_argument('-s', '--sample-name',
                            dest='sample_name',
                            required=True,
                            help='Sample name')
        parser.add_argument('-r', '--region',
                            dest='region',
                            default=None,
                            help='Region')
        parser.add_argument('-e', '--extension-bp',
                            dest='ext',
                            required=True,
                            type=int,
                            help='Extended bp for normal junctions')
        parser.add_argument('-p', '--ploidy',
                            dest='ploidy',
                            required=True,
                            default=2,
                            type=int,
                            help='Extended bp for normal junctions')
        parser.add_argument('-c', '--out-config',
                            dest='out_config',
                            required=True,
                            help='Output path of config')
        parser.add_argument('-g', '--segment',
                            dest='seg',
                            required=True,
                            help='Output path of segment')
        parser.add_argument('-i', '--keep-imprecise',
                            dest='keep_imprecise',
                            action='store_true',
                            default=True,
                            help='Keep imprecise SV')
        parser.add_argument('-I', '--keep-insertions',
                            dest='keep_insertions',
                            action='store_true',
                            default=True,
                            help='Keep insertions')
        args = parser.parse_args(sys.argv[2:])

        # chrom = args.region.split(':')[0]
        # start, end = [int(i) for i in args.region.split(':')[1].split('-')]
        bps_map = pd.read_csv(args.bps_map, sep='\t')
        chrom_infos = pd.read_csv(args.chrom_info, sep='\t')
        print('Reading SV')
        if not args.is_seeksv:
            sv = bpsmap.read_sv(args.sv_file)
            sv = bpsmap.get_precise_sv(sv, drop_imprecise=not args.keep_imprecise, drop_insertions=not args.keep_insertions)
    # sv = bpsmap.get_precise_sv_svaba(sv, chrom, start, end)
        # n, nc, sv = bpsmap.merge_sv_tgs2sgs(sv, sv, 10)
        else:
            sv = pd.read_table(args.sv_file, skiprows=1, header=None,
                                  usecols=[0, 1, 2, 3, 4, 5, 6, 7],
                                  names=['chrom_5p', 'pos_5p', 'strand_5p', 'left_read',
                                         'chrom_3p', 'pos_3p', 'strand_3p', 'right_read'])
        # config.map_bps_sv(sv, bps_map)
        config.map_bps_chrom_infos(chrom_infos, bps_map)
        sv = config.dedup(sv)

        segs = pd.DataFrame()
        id_start = 1
        for row in chrom_infos.itertuples():
            seg, id_start = config.segmentation(sv, row.chrom, row.start, row.end, id_start, drop_imprecise=not args.keep_imprecise, drop_insertions=not args.keep_insertions)
            segs = segs.append(seg)

        # segs = config.segmentation(sv, chrom, start, end)
        segs.to_csv(args.seg, index=False, sep='\t')
        bam = pysam.AlignmentFile(args.bam_file)
        depth_tabix = pysam.TabixFile(args.depth_file)
        # print('Calculating avg depth')
        # avg_depth = config.get_avg_depth(depth_tabix, chrom, start, end)
        # avg_depth = 300

        print('Updating junc db')
        # if args.junc_db:
        #     junc_db = pd.read_table(args.junc_db)
        #     junc_db = config.update_junc_db_by_sv(sv, junc_db, depth_tabix)
        #     junc_db = config.update_junc_db_by_seg(segs, junc_db, depth_tabix, chrom, start, end, bam, args.ext)
        #     config.write_junc_db(args.junc_db, junc_db)
        # if args.new_junc_db:
        junc_db = pd.DataFrame(columns=['chrom_5p', 'pos_5p', 'strand_5p', 'chrom_3p', 'pos_3p', 'strand_3p', 'count'])
        junc_db = config.update_junc_db_by_sv(sv, junc_db)
        # for chrom in segs.chrom.unique():
        #     chrom_info = chrom_infos.loc[lambda row: row.chrom == chrom].iloc[0]
        #     seg = segs.loc[lambda row: row.chrom == chrom]
        junc_db = config.update_junc_db_by_seg_in_chrom(segs, junc_db, bam, args.ext)
        config.write_junc_db(args.junc_db, junc_db)

        config.generate_config(args.out_config, args.sample_name, sv, segs, depth_tabix, bam, ext=args.ext, ploidy=args.ploidy)

    def mergedb(self):
        import pandas as pd
        import config
        parser = argparse.ArgumentParser(description='Merge junction DB')
        parser.add_argument('-l', '--junc-list',
                            dest='junc_list',
                            required=True,
                            help='List of junction DB')
        parser.add_argument('-J', '--seg-junc',
                            dest='seg_junc',
                            default=None,
                            help='Segment junction')
        parser.add_argument('-b', '--bps-map',
                            dest='bps_map',
                            required=True,
                            help='Break point map')
        parser.add_argument('-o', '--out-db',
                            dest='out_db',
                            required=True,
                            help='Output DB')
        args = parser.parse_args(sys.argv[2:])

        bps_map = pd.read_csv(args.bps_map, sep='\t')
        if args.seg_junc is None:
            seg_junc = None
        else:
            seg_junc = pd.read_csv(args.seg_junc, sep='\t',
                                     names=[
                                         'chrom_5p', 'pos_5p', 'strand_5p',
                                         'chrom_3p', 'pos_3p', 'strand_3p', 'count'
                                         ]
            )
            config.map_bps_junc(seg_junc, bps_map)
        merged = pd.DataFrame(columns=['chrom_5p', 'pos_5p', 'strand_5p', 'chrom_3p', 'pos_3p', 'strand_3p', 'count'])
        with open(args.junc_list, 'r') as fin:
            n = 1
            for line in fin:
                line = line[:-1]
                # print(f'{n} {line}')
                junc = pd.read_csv(line, sep='\t')
                for row in junc.itertuples():
                    idx1 = merged.loc[lambda r: r.chrom_5p == row.chrom_5p]\
                                 .loc[lambda r: r.pos_5p == row.pos_5p]\
                                 .loc[lambda r: r.strand_5p == row.strand_5p]\
                                 .loc[lambda r: r.chrom_3p == row.chrom_3p]\
                                 .loc[lambda r: r.pos_3p == row.pos_3p]\
                                 .loc[lambda r: r.strand_3p == row.strand_3p].index
                    if idx1.empty:
                        merged = merged.append({'chrom_5p': row.chrom_5p,
                                                'pos_5p': row.pos_5p,
                                                'strand_5p': row.strand_5p,
                                                'chrom_3p': row.chrom_3p,
                                                'pos_3p': row.pos_3p,
                                                'strand_3p': row.strand_3p,
                                                'count': row.count}, ignore_index=True)
                    else:
                        merged.at[idx1[0], 'count'] += row.count
                n += 1
        if seg_junc is not None:
            for row in seg_junc.itertuples():
                idx1 = merged.loc[lambda r: r.chrom_5p == row.chrom_5p]\
                             .loc[lambda r: r.pos_5p == row.pos_5p]\
                             .loc[lambda r: r.strand_5p == row.strand_5p]\
                             .loc[lambda r: r.chrom_3p == row.chrom_3p]\
                             .loc[lambda r: r.pos_3p == row.pos_3p]\
                             .loc[lambda r: r.strand_3p == row.strand_3p].index
                if idx1.empty:
                    merged = merged.append({'chrom_5p': row.chrom_5p,
                                            'pos_5p': row.pos_5p,
                                            'strand_5p': row.strand_5p,
                                            'chrom_3p': row.chrom_3p,
                                            'pos_3p': row.pos_3p,
                                            'strand_3p': row.strand_3p,
                                            'count': row.count}, ignore_index=True)
                else:
                    merged.at[idx1[0], 'count'] += row.count
        merged.sort_values(by=['chrom_5p', 'pos_5p', 'strand_5p', 'count']).to_csv(args.out_db, sep='\t', index=False)

    def parseILP(self):
        import parseILP
        parser = argparse.ArgumentParser(description='Parse ILP results')
        parser.add_argument('-i', '--in-checked-lh',
                            dest='checked_lh',
                            required=True,
                            help='Checked lh file')
        parser.add_argument('-s', '--sol',
                            dest='sol_file',
                            required=True,
                            help='Solution file')
        parser.add_argument('-o', '--output-balanced-lh',
                            dest='balanced_lh',
                            required=True,
                            help='Output balanced lh file')
        args = parser.parse_args(sys.argv[2:])

        sol = parseILP.parse_ilp_result(args.sol_file)
        parseILP.generate_balanced_lh(args.balanced_lh, args.checked_lh, sol)

if __name__ == '__main__':
    MainArgParser()
