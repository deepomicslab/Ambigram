from pyfaidx import Fasta,Faidx
import sys
import os
import random
import argparse
import logging
import e_size
import re

PYTHON = "~/miniconda3/envs/py3/bin/python"
LOCALHAP = "/home/gzpan2/app/localhaptgs/debug/localHap"
SAMTOOLS = "~/app/samtools/bin/samtools"
CBC = "~/miniconda3/envs/py2/bin/cbc"
pbsim = "~/app/pbsim2/src/pbsim"
pbmodel = "~/app/pbsim2/data/P6C4.model"
hpvpip_root = "~/app/hpvpip"
faToTwoBit = "~/app/faToTwoBit"
computeGCBias = "/home/xuedowang2/app/conda/envs/py37/bin/computeGCBias"
correctGCBias = "/home/xuedowang2/app/conda/envs/py37/bin/correctGCBias"
samtools = "~/app/samtools/bin/samtools"
# ~/app/pbsim2/src/pbsim --depth 20 --prefix tgs --hmm_model ~/app/pbsim2/data/P6C4.model test.out.fa
def execmd(cmd):
    print("Exec: {}".format(cmd))
    logging.info(cmd)
    # os.system(cmd)
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--host_ref',
                        required=True,
                        help='host_ref')
    parser.add_argument('--v_ref',
                        required=True,
                        help='v_ref')
    parser.add_argument('--out',
                        required=True,
                        help='out')
    parser.add_argument('--script_root',
                        required=True,
                        help='out')
    parser.add_argument('--mutforge',
                        required=True,
                        help='mutforge')
    parser.add_argument('--simple_par',
                        required=True,
                        help='simple_par')
    parser.add_argument('--depth',
                        required=True,
                        help='simple_par')
    parser.add_argument('--host_count',
                        required=True,
                        type=int,
                        help='simple_par')
    parser.add_argument('--s_times',
                        required=True,
                        type=int,
                        help='simple_par')
    args = parser.parse_args()
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    # choose host chr
    # mk_fa(args.host_ref, host_chrs, args.v_ref, args.v_chr, args.out)
    # generate_var(host_chrs, args.v_chr, v_len,args.out,args.v_ref)
    simulate(args.mutforge, args.host_ref, args.v_ref,args.simple_par, args.out, args.script_root, args.depth, args.host_count, args.s_times)
# def sim_tgs(out_dir,depth = 20):
#     cmd1 = "{} --depth {} --prefix tgs --hmm_model {} test.out.fa"

def g_tgs_ref(out_dir,all_chrs, depth = 20):
    out_fa_file = out_dir+".out.fa"
    out_fa = open(out_fa_file,"w")
    res = {}
    for i in all_chrs:
        res[i] = []
    ref_host = Fasta(out_dir+".hap1.fa")
    for k in list(ref_host.keys()):
        if "original" in k:
            continue
        tmp_chr = k.split(":")[5].split("=")[1]
        res[tmp_chr].append(k)
    for key,v in res.items():
        combined_fa = ""
        for k in v:
            combined_fa = combined_fa+str(ref_host[k][0:])
        out_fa.write(">{}\n".format(key))
        out_fa.write(combined_fa+"\n")
    out_fa.close()
#     simulate tgs
    cmd1 = "{} --depth {} --prefix {} --hmm_model {} {}".format(pbsim,depth, out_dir,pbmodel, out_fa_file)
    cmd2 = "cat {}_*.fastq > {}.tgs.fastq".format(out_dir,out_dir)
    cmd3 = "sed -n '1~4s/^@/>/p;2~4p' {}.tgs.fastq > {}.tgs.fasta".format(out_dir, out_dir)
    cmd4 = "{} {}/main.py process_tgs --ref {}/mix.fa -l {}.lh -t {}.tgs.fasta -o {}/tgsout --max_bias 0.2".format(PYTHON, hpvpip_root, out_dir,out_dir,out_dir,out_dir)
    execmd(cmd1)
    execmd(cmd2)
    execmd(cmd3)
    execmd(cmd4)

def parse_mean_depth(bam,out_dir,n_size):
    cmd = "{} coverage {} > {}.scov".format(samtools, bam,out_dir)
    execmd(cmd)
    for l in open(out_dir+".scov") :
        if "#" not in l:
            a = re.split("\s+",l)
            o_size = int(a[2])
            o_depth = int(a[6])
            n_depth = o_depth*(o_size/n_size)
            break
    return n_depth
def check_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir,exist_ok=True)

def gc_correction(input_bam, out_dir, effectiveGenomeSize):
    ref = out_dir + "/mix.fa"
    corrected_bam = out_dir+".gc.bam"
    # faToTwoBit hg38_hpv.fa hg38_hpv.bit
    cmd1 = "echo 'skip ref faToTwoBit'"
    cmd2 = "echo 'skip generate freq.txt'"
    if not check_dir(ref+".freq.txt"):
        cmd1 = "{} {} {}.2bit".format(faToTwoBit, ref, ref)
        cmd2 = "{} -b {} --effectiveGenomeSize {} -g {}.2bit --GCbiasFrequenciesFile {}.freq.txt".format(computeGCBias, input_bam, effectiveGenomeSize, ref, ref)
    cmd3 = "{} -b {} --effectiveGenomeSize {} -g {}.2bit --GCbiasFrequenciesFile {}.freq.txt -o {}".format(correctGCBias, input_bam, effectiveGenomeSize, ref, ref, corrected_bam)
    cmd4 = "{} index {} -@ {}".format(samtools, corrected_bam, 48)
    execmd(cmd1)
    execmd(cmd2)
    execmd(cmd3)
    execmd(cmd4)
    return corrected_bam

def run_local(out_dir,script_root,vc,v_len,selected_chrs,depth, gc_bam):
    cmd_seek="bash {}/seek.sh {} {}.lib1.bam {}/mix.fa".format(script_root,out_dir,out_dir,out_dir)
    cmd_bps = "{} {}/main.py bpsmap -l {}.seek.sv.txt -o {} -v {} --v_len {} --h_chrs {} --out_bed {}.bed".format(PYTHON, script_root,out_dir,out_dir,vc,v_len, ','.join(selected_chrs),out_dir)
    bed_file = out_dir+".bed"
    cmd_depth= "{} depth -aa -b {}.bed {}.gc.bam | bgzip -c > {}.depth.gz && tabix -s 1 -b 2 -e 2 {}.depth.gz".format(SAMTOOLS,out_dir,out_dir,out_dir,out_dir)
    cmd_config = "{} {}/main.py config -f {}.seek.sv.txt -b {} -m {}.bps -j {}.junc -d {}.depth.gz  -s {} -e 5 -c {}.lh -g {}.seg -p 2 -S --v_chr {} --avg_whole_dp {}\
         --v_len {} --h_chrs {}".format(PYTHON,script_root, out_dir,gc_bam,out_dir,out_dir,out_dir,out_dir,out_dir,out_dir,vc,depth,v_len,','.join(selected_chrs))
    cmd_check = "{} --op check --juncdb {}.junc --in_lh {}.lh --out_lh {}.checked.lh --lp_prefix {} --verbose".format(LOCALHAP, out_dir, out_dir, out_dir, out_dir)
    cmd_cbc = "{} {}.lp solve solu {}.sol".format(CBC,out_dir,out_dir)
    cmd_parse = "{} {}/main.py parseILP -i {}.checked.lh -s {}.sol -o {}.balanced.lh".format(PYTHON, script_root,out_dir,out_dir,out_dir)
    cmd_solve = "{} --op solve --juncdb {}.junc --in_lh {}.balanced.lh --circuits {}.circuits --hap {}.haps --traversed {}.traversed --tgs_order {}/tgsout/tgs.juncs --verbose".format(LOCALHAP,out_dir,out_dir,out_dir,out_dir,out_dir,out_dir)
    execmd(cmd_seek)
    execmd(cmd_bps)
    execmd(cmd_depth)
    execmd(cmd_config)
    execmd(cmd_check)
    execmd(cmd_cbc)
    execmd(cmd_parse)
    selected_chrs.append(vc)
    g_tgs_ref(out_dir,selected_chrs)
    execmd(cmd_solve)


def simulate(mutforge, host_ref, v_ref,simple_par, out, script_root,depth, host_count, s_times):
    ref_host = Fasta(host_ref)
    ref_v = Fasta(v_ref)
    host_chrs = list(ref_host.keys())[0:-3]
    v_chrs = {}
    total_size = 0
    for k in list(ref_v.keys()):
        v_chrs[k] = ref_v[k][0:].end
    for vc,v_len in v_chrs.items():
        for i in range(0,s_times):
            total_size = 0
            # selected_chrs = random.choices(host_chrs, k=3)
            # if i > 2:
            selected_chrs = random.choices(host_chrs, k=host_count)
            while len(set(selected_chrs)) != host_count:
                selected_chrs = random.choices(host_chrs, k=host_count)
            # if i >= 5:
            #     selected_chrs = random.choices(host_chrs, k=3)
            #     while len(set(selected_chrs)) != 3:
            #         selected_chrs = random.choices(host_chrs, k=3)
            out_dir = os.path.join(out, "{}_{}".format(vc,"_".join(selected_chrs)),"test")
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            mix_fa = mk_fa(host_ref, selected_chrs, v_ref, vc, out_dir)
            var_file = generate_var(selected_chrs, vc, v_len, out_dir, v_ref)
            cmd_mu = "{} -x bam -n 128 -v {} {} {} {} -o {} -b 0".format(mutforge, var_file, mix_fa, simple_par,mix_fa,out_dir)
            # cmd_seek="bash seek.sh {} {}.lib1.bam {}/mix.fa".format(out_dir,out_dir,out_dir)
            # cmd_bps = "{} {}/main.py bpsmap -l {}.seek.sv.txt -o {}.bps -v {} --v_len {} --h_chrs {}".format(python, script_root,out_dir,out_dir,vc,v_len, ','.join(selected_chrs))
            # cmd_depth= "{} depth -b {} {}.lib1.bam > {}.txt".format(samtools,bed_file,out_dir,out_dir)
            # cmd_config = "{} "
            # cmd_check = "{} --op check --juncdb {}.junc --in_lh {}.lh --out_lh {}.checked.lh --lp_prefix {} --verbose".format(localhap, out_dir, out_dir, out_dir, out_dir)
            # cmd_cbc = "{} {}.lp solve solu {}.sol".format(cbc,out_dir,out_dir)
            # cmd_parse = "{} {}/main.py parseILP -i {}.checked.lh -s {}.sol -o {}.balanced.lh".format(python, script_root,out_dir,out_dir,out_dir)
            # cmd_solve = "{} --op solve --juncdb {}.junc --in_lh {}.balanced.lh --circuits {}.circuits --hap {}.haps --verbose".format(localhap,out_dir,out_dir,out_dir,out_dir)
            execmd(cmd_mu)
            total_size = total_size + v_len + sum([e_size.sizes[c] for c in selected_chrs])
            gc_bam = gc_correction(out_dir+".lib1.bam",out_dir,total_size)
            n_depth = parse_mean_depth(gc_bam,out_dir,total_size)
            run_local(out_dir,script_root,vc,v_len,selected_chrs,n_depth,gc_bam)
            # execmd(cmd_seek)
def mk_fa(host_ref,host_chrs,v_ref,v_chr,out):
    # extract
    ref_host = Fasta(host_ref)
    ref_v = Fasta(v_ref)
    mix_f = os.path.join(out,"mix.fa")
    mix_out = open(mix_f,"w")
    for el in host_chrs:
        mix_out.write(">{}\n".format(el))
        mix_out.write(str(ref_host[el][0:])+"\n")
    mix_out.write(">{}\n".format(v_chr))
    mix_out.write(str(ref_v[v_chr][0:]))
    mix_out.close()
    fa = Faidx(mix_f)
    # cmd = "bwa index {}".format(mix_f)
    # execmd(cmd)
    return mix_f
def in_region(i,rs):
    for r in rs:
        if r[0]<=i and i <= r[1]:
            return True;
    return False

def generate_var(host_chrs,v_chr,v_len,out,fa_file):
    var_file=os.path.join(out,"mix.var")
    var_out = open(var_file, "w")
    chr_num = len(host_chrs)
    for hc in host_chrs:
        r_start = random.randint(20000000, 25005000)
        r_end = random.randint(r_start + 1000, r_start+10000)
        del_regions = []
        for i in range(0,3):
            tmp0 = ""
            tmp1 = ""
            pos = random.randint(r_start, r_end)
            # del
            if pos%3==0:
                d_size = random.randint(500,3000)
                while in_region(pos, del_regions) or in_region(pos+d_size,del_regions):
                    pos = random.randint(r_start, r_end)
                    d_size = random.randint(500,3000)
                del_regions.append([pos, pos+d_size])
                tmp0 = "VAR_{}_0_{}\tFDEL_{}\t1\t0\t{}\t{}\tTrue\t{}\tFalse\tNone".format(hc,i,i,hc,pos,d_size)
                tmp1 = "VAR_{}_1_{}\tFDEL_{}\t1\t1\t{}\t{}\tTrue\t{}\tFalse\tNone".format(hc,i,i,hc,pos,d_size)
            else:
                # ins
                rev = "f"
                i_start = random.randint(700,v_len - 500)
                i_end = random.randint(i_start + 300,i_start + 700)
                while in_region(i_start, del_regions) or in_region(i_end,del_regions):
                    i_start = random.randint(100,v_len - 1500)
                    i_end = random.randint(i_start + 500,v_len)
                times = random.randint(1,3)
                if i_start%3 == 0:
                    rev = "r"
                tmp0 = "VAR_{}_0_{}\tFINS_{}\t1\t0\t{}\t{}\tFalse\t0\tTrue\t{},{}:{}-{},{},{}".format(hc,i,i,hc,pos,fa_file,v_chr,i_start,i_end,times,rev)
                tmp1 = "VAR_{}_1_{}\tFINS_{}\t1\t1\t{}\t{}\tFalse\t0\tTrue\t{},{}:{}-{},{},{}".format(hc,i,i,hc,pos,fa_file,v_chr,i_start,i_end,times,rev)
            var_out.write(tmp0+"\n")
            var_out.write(tmp1+"\n")
    return var_file

if __name__ == "__main__":
    main()
