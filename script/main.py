import argparse

import os
from re import S
import sys

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
        d = d.append(sv.query(
            f'chrom_3p.str.contains("{args.ref_common_str}") & chrom_5p.str.contains("{args.ref_common_str}")'))
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
        parser.add_argument('-v', '--v_chr',
                            dest='v_chr',
                            default=None,
                            help='Region')
        parser.add_argument('--h_chrs',
                            dest='h_chrs',
                            default=None,
                            help='Region')
        parser.add_argument('--v_len',
                            dest='v_len',
                            default=None,
                            help='Region')
        parser.add_argument('-S', '--chrom-sv',
                            dest='chrom_sv',
                            required=False,
                            default=None,
                            help='Chromosome SV')
        # parser.add_argument('-C', '--chrom-info',
        #                     dest='chrom_info',
        #                     required=True,
        #                     help='Chromosome information')
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
        parser.add_argument('--out_bed',
                            default=True,
                            help='Keep insertions')
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
        seek_sv = True
        if seek_sv:
            sv_df = pd.read_table(args.sv_list, skiprows=1, header=None,
                                  usecols=[0, 1, 2, 3, 4, 5, 6, 7],
                                  names=['chrom_5p', 'pos_5p', 'strand_5p', 'left_read',
                                         'chrom_3p', 'pos_3p', 'strand_3p', 'right_read'])
            print(sv_df['pos_3p'])
            sv_df = sv_df.astype({
                'chrom_5p':str, 'pos_5p':np.int64, 'strand_5p':str, 'left_read':str,
                'chrom_3p':str, 'pos_3p':np.int64, 'strand_3p':str, 'right_read':str,
            })
            # sv_df = pd.read_csv(args.sv_list)
        else:
            # sv_df = bpsmap.concat_sv(args.sv_list)
            sv_df = bpsmap.parse_sur(args.sv_list)
        chroms = args.h_chrs.split(",")
        chroms.append(args.v_chr)
        # chrom_infos = pd.read_csv(args.chrom_info, dtype={'end': np.int64}, sep='\t')
        # print(chrom_junc.head())
        # print(sv_df.head())
        # print(chrom_infos.head())
        if args.region:
            ## TODO: if region not None
            pass
        else:
            chroms = chroms
            if chrom_junc is not None:
                chroms = sorted(set(
                    list(chrom_junc.chrom_5p.unique()) +
                    list(chrom_junc.chrom_3p.unique()) +
                    chroms))
            else:
                # print(chroms)
                print(list(sv_df.columns))
            print(chroms,"xxx")
            chroms = sorted(set(chroms))
            bps_map = []
            print(sv_df['chrom_5p'])
            for chrom in chroms:
                sv_5p = bpsmap.get_precise_sv(sv_df, chrom_5p=chrom, drop_imprecise=not args.keep_imprecise,
                                              drop_insertions=not args.keep_insertions)
                sv_3p = bpsmap.get_precise_sv(sv_df, chrom_3p=chrom, drop_imprecise=not args.keep_imprecise,
                                              drop_insertions=not args.keep_insertions)
                if chrom_junc is not None:
                    chrom_junc_5p = chrom_junc.loc[lambda df: df.chrom_5p == chrom]
                    chrom_junc_3p = chrom_junc.loc[lambda df: df.chrom_3p == chrom]
                print(chrom)
                print((sv_3p))
                bps = ""
                if chrom == args.v_chr:
                    bps = bpsmap.get_breakpoints(sv_5p, sv_3p, True) + [1, int(args.v_len)]
                else:
                    bps = bpsmap.get_breakpoints(sv_5p, sv_3p, False)
                print(bps)
                if chrom_junc is not None:
                    bps = np.array(sorted(set(
                        bps +
                        bpsmap.get_breakpoints(chrom_junc_5p, chrom_junc_3p)
                    )))
                else:
                    bps = np.array(sorted(set(bps)))
                bps_map.extend([(chrom, *t) for t in bpsmap.map_bps(bps, 10)])
            bpsmap.write_bps_map(args.bps_map_out, bps_map,chroms)
            # bpsmap.generate_bed(bps_map,args.out_bed)

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
        parser = argparse.ArgumentParser(description='Generate localhap config for each individual')
        parser.add_argument('-f', '--sv-file',
                            dest='sv_file',
                            required=True,
                            help='Individual SV file')
        parser.add_argument('-b', '--bam-file',
                            dest='bam_file',
                            required=True,
                            help='Individual BAM file')
        parser.add_argument('--v_chr',
                            dest='v_chr',
                            required=True,
                            help='Individual BAM file')
        parser.add_argument('--h_chrs',
                            dest='h_chrs',
                            required=True,
                            help='Individual BAM file')
        parser.add_argument('--v_len',
                            dest='v_len',
                            type=int,
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
        parser.add_argument('--avg_whole_dp',
                            dest='avg_whole_dp',
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
        chroms = args.h_chrs.split(",")
        chroms.append(args.v_chr)

        # chrom = args.region.split(':')[0]
        # start, end = [int(i) for i in args.region.split(':')[1].split('-')]
        bps_map = pd.read_csv(args.bps_map, sep='\t')
        print('Reading SV')
        if not args.is_seeksv:
            # sv = bpsmap.read_sv(args.sv_file)
            sv = bpsmap.parse_sur(args.sv_file)
            sv = bpsmap.get_precise_sv(sv, drop_imprecise=not args.keep_imprecise,
                                       drop_insertions=not args.keep_insertions)
        # sv = bpsmap.get_precise_sv_svaba(sv, chrom, start, end)
        # n, nc, sv = bpsmap.merge_sv_tgs2sgs(sv, sv, 10)
        else:
            sv = pd.read_table(args.sv_file, skiprows=1, header=None,
                               usecols=[0, 1, 2, 3, 4, 5, 6, 7],
                               names=['chrom_5p', 'pos_5p', 'strand_5p', 'left_read',
                                      'chrom_3p', 'pos_3p', 'strand_3p', 'right_read'])
        config.map_bps_sv(sv, bps_map)
        # config.map_bps_chrom_infos(chrom_infos, bps_map)
        sv = config.dedup(sv)

        segs = pd.DataFrame()
        beds = []
        id_start = 1
        for chrom in chroms:
            seg, id_start,i_start,i_end = config.segmentation(sv, chrom,args.v_chr,args.v_len, id_start,
                                                drop_imprecise=not args.keep_imprecise,
                                                drop_insertions=not args.keep_insertions)
            segs = segs.append(seg)
            beds.append([chrom, str(i_start),str(i_end)])

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

        config.generate_config(args.out_config, args.sample_name, sv, segs, depth_tabix, bam,args.v_chr,args.avg_whole_dp, ext=args.ext,
                               ploidy=args.ploidy)

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
                    idx1 = merged.loc[lambda r: r.chrom_5p == row.chrom_5p] \
                        .loc[lambda r: r.pos_5p == row.pos_5p] \
                        .loc[lambda r: r.strand_5p == row.strand_5p] \
                        .loc[lambda r: r.chrom_3p == row.chrom_3p] \
                        .loc[lambda r: r.pos_3p == row.pos_3p] \
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
                idx1 = merged.loc[lambda r: r.chrom_5p == row.chrom_5p] \
                    .loc[lambda r: r.pos_5p == row.pos_5p] \
                    .loc[lambda r: r.strand_5p == row.strand_5p] \
                    .loc[lambda r: r.chrom_3p == row.chrom_3p] \
                    .loc[lambda r: r.pos_3p == row.pos_3p] \
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
    
    # change sv positions in bed.txt
    def updateBed(self):
        parser = argparse.ArgumentParser(description='Change sv positions in bed.txt according to sv file')
        parser.add_argument('-i', '--in_sv',
                            dest='sv',
                            required=True,
                            help='Input sv file')
        parser.add_argument('-b', '--bed_file',
                            dest='bed',
                            required=True,
                            help='A bed file constructed by localHap bfb path.')
        args = parser.parse_args(sys.argv[2:])
        #read sv infomation
        sv_file = open(args.sv, "r")
        sv = []
        for line in sv_file.readlines():
            if line[:8] == "chrom_5p":
                continue
            arr = line.split("\t")
            arr[1] = int(arr[1])
            arr[4] = int(arr[4])
            sv.append(arr)
            print(sv[-1])
        #read bed strings
        bed_file = open(args.bed, "r")
        bed = []
        for line in bed_file.readlines():
            arr = line.split(" ")
            arr[1] = int(arr[1])
            arr[2] = int(arr[2])
            arr[-1] = arr[-1][0]
            bed.append(arr)
            print(bed[-1])
        #update break positions
        for i in range(0, len(bed)-1):
            if bed[i][0] == bed[i+1][0] and bed[i][3] == bed[i+1][3]:
                continue
            #fold-back inversion
            for info in sv:
                if bed[i][0] in info and bed[i+1][0] in info:
                    if  (info[1] in range(bed[i][1], bed[i][2]+1) and info[4] in range(bed[i+1][1], bed[i+1][2]+1)) or \
                        (info[4] in range(bed[i][1], bed[i][2]+1) and info[1] in range(bed[i+1][1], bed[i+1][2]+1)): 
                        pos1 = 0
                        pos2 = 0
                        if info[2] == bed[i][-1] and info[5] == bed[i+1][-1]:
                            pos1 = info[1]
                            pos2 = info[4]
                        elif info[5] == bed[i][-1] and info[2] == bed[i+1][-1]:
                            pos1 = info[4]
                            pos2 = info[1]
                        else:
                            continue
                        #for translocation
                        if info[0] != info[3]:
                            if info[1] in range(bed[i][1], bed[i][2]+1):
                                pos1 = info[1]
                                pos2 = info[4]
                            else:
                                pos1 = info[4]
                                pos2 = info[1]
                        #change the positions
                        if bed[i][3] == "forward":
                            bed[i][2] = pos1
                        else:
                            bed[i][1] = pos1
                        if bed[i+1][3] == "forward":
                            bed[i+1][1] = pos2
                        else:
                            bed[i+1][2] = pos2
        #re-write the bed file
        res = ""
        for info in bed:
            for entry in info:
                res += str(entry)+" "
            res += "\n"
        bed_file = open(args.bed, "w")
        bed_file.write(res)  

    def bfb2fasta(self):
        import pybedtools
        parser = argparse.ArgumentParser(description='Convert a segment (in bed format) into a bp sequence according to the fasta file')
        parser.add_argument('-i', '--in_fasta',
                            dest='inFasta',
                            required=True,
                            help='Input fasta file')
        parser.add_argument('-b', '--bed_file',
                            dest='bed',
                            required=True,
                            help='A bed file constructed by bfb path.')
        # parser.add_argument('-s', '--segment',
        #                     dest='segment',
        #                     required=True,
        #                     help='A string of segment in bed format. e.g., chr6:150336:256499')
        parser.add_argument('-t', '--out_txt',
                            dest='outTxt',
                            required=True,
                            help='A string of output directory')
        parser.add_argument('-o', '--out_fasta',
                            dest='outFasta',
                            required=True,
                            help='A fasta file consisting of bfb path')
        args = parser.parse_args(sys.argv[2:])
        #read bed strings
        bed_file = open(args.bed, "r")
        lines = bed_file.readlines()
        bedStr = ""
        for line in lines:
            bedStr += line+"\n"
        #extract the bp sequence from fasta
        seg = pybedtools.BedTool(bedStr, from_string=True)
        fasta = pybedtools.example_filename(args.inFasta)
        seg = seg.sequence(fi=fasta, s=True)
        #output the result
        # if os.path.exists(args.output):
        txt = open(args.outTxt, "w")
        txt.write(open(seg.seqfn).read())
        # print(open(seg.seqfn).read())
        #convert BedTool output into bp sequence (with inversion and translocation)
        txt = open(args.outTxt, "r")
        lines = txt.readlines()
        res = ">BFBPATH\n"
        for i in range(0, len(lines)):
            if i%2 != 0:
                continue
            res += lines[i+1][:-1]#the last char: \n
            print(i/2+1, len(lines[i+1]))
            # print(lines[i+1][:-1])

        outFile = open(args.outFasta, "w")
        outFile.write(res)

    def vcf2sv(self):
        parser = argparse.ArgumentParser(description='Convert a vcf.txt file into a sv file for generate_lh()')
        parser.add_argument('-i', '--in_vcf',
                            dest='vcf',
                            required=True,
                            help='Input vcf file')
        parser.add_argument('-o', '--out_sv',
                            dest='output',
                            required=True,
                            help='A name of output')
        parser.add_argument('-b', '--barcode_group',
                            dest='barcode',
                            required=False,
                            help='A barcode_group.csv file')
        parser.add_argument('-g', '--group',
                            dest='group',
                            required=False,
                            help='A group_meta.csv file')
        args = parser.parse_args(sys.argv[2:])
        vcf = open(args.vcf, "r")
        if args.barcode == None:            
            sv, arr = [], []
            for line in vcf.readlines():            
                entry = line.split("\t")
                # if entry[0]==entry[3] and (entry[2]==entry[5] or abs(int(entry[1])-int(entry[4])) > 100000):
                #     continue
                # if int(entry[7])+int(entry[8])<20 or int(entry[8])<6:
                #     continue                
                depth = entry[13].split('DP:')[1].split(' ')[0]
                arr.append([entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], depth])
            # filter noisy SV
            if arr[1][0] == arr[1][3] and arr[1][2] != arr[1][5]:
                sv.append(arr[0])
            for i in range(1, len(arr)-1):
                if arr[i][0] != arr[i][3] or arr[i][2] != arr[i][5]:
                    sv.append(arr[i])
                else:
                    if arr[i-1][0] == arr[i-1][3] and arr[i-1][2] != arr[i-1][5] or \
                        arr[i+1][0] == arr[i+1][3] and arr[i+1][2] != arr[i+1][5]:
                        sv.append(arr[i])
            if arr[-2][0] == arr[-2][3] and arr[-2][2] != arr[-2][5]:
                sv.append(arr[-1])
            # sv = sorted(sv, key=lambda d: int(d[1]))
            res = "chrom_5p\tbkpos_5p\tstrand_5p\tchrom_3p\tbkpos_3p\tstrand_3p\tavg_cn\n"
            for d in sv:
                res += '\t'.join(d)+'\n'
            outFile = open('{}_sv.txt'.format(args.output), "w")
            outFile.write(res)
        else:
            barcode, sv = [], []
            for line in vcf.readlines():            
                entry = line.split("\t")
                # if entry[2]==entry[5] or float(entry[10]) < 20:
                #     continue
                # if entry[0]==entry[3] and entry[2]==entry[5]:
                #     continue
                codes = []
                for code in entry[13].split(';')[0].split(','):
                    if len(codes) == 0:
                        codes.append(code[3:-2])
                    else:
                        codes.append(code[:-2])
                barcode.append(codes)
                str1, str2 = 'h', 't'
                if entry[2] == '-':
                    str1 = 't'
                if entry[5] == '-':
                    str2 = 'h'
                depth = entry[13].split('DP:')[1].split(' ')[0]
                sv_str = entry[0]+"\t"+entry[1]+"\t"+entry[2]+"\t"
                sv_str += entry[3]+"\t"+entry[4]+"\t"+entry[5]+"\t"
                sv_str += depth #+"\t"+(str1+str2)+'\tNone\tNone\t'
                sv.append(sv_str)
            # read barcode file
            barcode2group = {}
            for line in open(args.barcode, "r").readlines():
                if line[0] == ',':
                    continue
                info = line.strip('\n').split(',')
                if args.group == None:
                    barcode2group[info[0]] = info[1]
                else:
                    barcode2group[info[1]] = info[2]
            # read group file
            group2subclone = {}
            if args.group != None:
                for line in open(args.group, 'r').readlines():
                    info = line.strip('\n').split(',')
                    if info[0] == 'label':
                        continue
                    group2subclone[info[1]] = info[0]
            # construct subclone
            subclone, cloneCount = {}, [{} for i in range(len(barcode))]
            for i in range(len(barcode)):
                for code in barcode[i]:
                    if code not in barcode2group.keys():
                        continue
                    if args.group == None:
                        k = barcode2group[code]
                    else:
                        k = group2subclone[barcode2group[code]]
                    if k not in subclone.keys():
                        subclone[k] = []
                    subclone[k].append(i)
                    if k not in cloneCount[i].keys():
                        cloneCount[i][k] = 0
                    cloneCount[i][k] += 1
            # output subclone sv
            for key, val in subclone.items():
                res = 'chrom_5p\tbkpos_5p\tstrand_5p\tchrom_3p\tbkpos_3p\tstrand_3p\tavg_cn\n'
                for i in list(set(val)):
                    clone = ''
                    for k, v in cloneCount[i].items():
                        clone += k+'='+str(v)+','
                    res += sv[i]+'\n'#+clone[:-1]+'\n'
                outFile = open('{}_{}_sv.txt'.format(args.output, key), "w")
                outFile.write(res)
        

    def seg2fasta(self):
        import pybedtools
        parser = argparse.ArgumentParser(description='Convert a segment file into the fasta file')
        parser.add_argument('-i', '--in_seg',
                            dest='seg',
                            required=True,
                            help='Input segment file')
        parser.add_argument('-r', '--ref_fasta',
                            dest='ref',
                            required=True,
                            help='Reference fasta file')
        parser.add_argument('-o', '--out_fasta',
                            dest='fasta',
                            required=True,
                            help='A name/directory of output')
        args = parser.parse_args(sys.argv[2:])
        #read segments
        segs = open(args.seg, "r")
        bedStr = ""
        for line in segs.readlines():
            info = line.split("\t")[0]
            chr = info.split(":")[0]
            pos = info.split(":")[1]
            bedStr += chr+" "+pos.split("-")[0]+" "+pos.split("-")[1]+" forward 1 +\n"
        print(bedStr)
        #extract the bp sequence from fasta
        seg = pybedtools.BedTool(bedStr, from_string=True)
        fasta = pybedtools.example_filename(args.ref)
        seg = seg.sequence(fi=fasta, s=True)
        #output the result
        outFasta = open(args.fasta, "w")
        outFasta.write(open(seg.seqfn).read())

    def parse_snif_vcf(self):
        parser = argparse.ArgumentParser(description='Convert a Sniffles vcf file into sv file')
        parser.add_argument('-i', '--in_vcf',
                            dest='vcf',
                            required=True,
                            help='Input vcf file')
        parser.add_argument('-o', '--out_sv',
                            dest='sv',
                            required=True,
                            help='Output sv file')
        args = parser.parse_args(sys.argv[2:])
        # read vcf file
        vcf = open(args.vcf, "r")
        res = "chrom_5p\tbkpos_5p\tstrand_5p\tchrom_3p\tbkpos_3p\tstrand_3p\tavg_cn\n"
        inv = ["++","--"]        
        trans = ["[","]","N"]
        for line in vcf.readlines():
            entry = line.split("\t")
            if line[0] != '#':
                strands = ""
                for item in entry[7].split(";"):
                    if item[:8] == "STRANDS=":
                        strands = item[-2:]
                        break                                             
                if strands in inv or (entry[4][0] in trans and len(entry[4])>1) or entry[4]=="<TRA>":
                    new_str = strands
                    if strands == "+-":
                        new_str = "++"
                    elif strands == "++":
                        new_str = "+-"
                    elif strands == "--":
                        new_str = "-+"
                    elif strands == "-+":
                        new_str = "--"
                    chrom_3p = entry[7].split(";")[2][5:]
                    pos_3p = entry[7].split(";")[3][4:]
                    num_vReads = entry[9].split(":")[-1]
                    if chrom_3p[0] != 'c':
                        chrom_3p = "chr"+chrom_3p
                    if entry[0][0] != 'c':
                        entry[0] = "chr"+entry[0]
                    res += entry[0]+"\t"+entry[1]+"\t"+new_str[0]+"\t" \
                        +chrom_3p+"\t"+pos_3p+"\t"+new_str[1]+"\t"+num_vReads
        #output the result
        bfb_sv = open(args.sv, "w")
        bfb_sv.write(res)

if __name__ == '__main__':
    MainArgParser()
