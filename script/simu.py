from pyfaidx import Fasta,Faidx
import sys
import os
import random
import argparse
PYTHON = "~/miniconda3/envs/py3/bin/python"
LOCALHAP = "/home/gzpan2/app/localhaptgs/debug"
SAMTOOLS = "~/app/samtools/bin/samtools"
CBC = "~/miniconda3/envs/py2/bin/cbc"
def execmd(cmd):
    print("Exec: {}".format(cmd))
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
    args = parser.parse_args()
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    # choose host chr
    # mk_fa(args.host_ref, host_chrs, args.v_ref, args.v_chr, args.out)
    # generate_var(host_chrs, args.v_chr, v_len,args.out,args.v_ref)
    simulate(args.mutforge, args.host_ref, args.v_ref,args.simple_par, args.out, args.script_root)

def run_local(out_dir,script_root,vc,v_len,selected_chrs):
    cmd_seek="bash {}/seek.sh {} {}.lib1.bam {}/mix.fa".format(script_root,out_dir,out_dir,out_dir)
    cmd_bps = "{} {}/main.py bpsmap -l {}.seek.sv.txt -o {} -v {} --v_len {} --h_chrs {} --out_bed {}.bed".format(PYTHON, script_root,out_dir,out_dir,vc,v_len, ','.join(selected_chrs),out_dir)
    bed_file = out_dir+".bed"
    cmd_depth= "{} depth -aa -b {}.bed {}.lib1.bam | bgzip -c > {}.depth.gz && tabix -s 1 -b 2 -e 2 {}.depth.gz".format(SAMTOOLS,out_dir,out_dir,out_dir,out_dir)
    cmd_config = "{} {}/main.py config -f {}.seek.sv.txt -b {}.lib1.bam -m {}.bps -j {}.junc -d {}.txt.gz  -s {} -e 5 -c {}.lh -g {}.seg -p 2 -S --v_chr {} --avg_whole_dp {}\
         --v_len {} --h_chrs {}".format(PYTHON,script_root, out_dir,out_dir,out_dir,out_dir,out_dir,out_dir,out_dir,out_dir,vc,40,v_len,','.join(selected_chrs))
    cmd_check = "{} --op check --juncdb {}.junc --in_lh {}.lh --out_lh {}.checked.lh --lp_prefix {} --verbose".format(LOCALHAP, out_dir, out_dir, out_dir, out_dir)
    cmd_cbc = "{} {}.lp solve solu {}.sol".format(CBC,out_dir,out_dir)
    cmd_parse = "{} {}/main.py parseILP -i {}.checked.lh -s {}.sol -o {}.balanced.lh".format(PYTHON, script_root,out_dir,out_dir,out_dir)
    cmd_solve = "{} --op solve --juncdb {}.junc --in_lh {}.balanced.lh --circuits {}.circuits --hap {}.haps --verbose".format(LOCALHAP,out_dir,out_dir,out_dir,out_dir)
    execmd(cmd_seek)
    execmd(cmd_bps)
    execmd(cmd_depth)
    execmd(cmd_config)
    execmd(cmd_check)
    execmd(cmd_cbc)
    execmd(cmd_parse)
    execmd(cmd_solve)


def simulate(mutforge, host_ref, v_ref,simple_par, out, script_root):
    ref_host = Fasta(host_ref)
    ref_v = Fasta(v_ref)
    host_chrs = list(ref_host.keys())[0:-3]
    v_chrs = {}
    for k in list(ref_v.keys()):
        v_chrs[k] = ref_v[k][0:].end
    for vc,v_len in v_chrs.items():
        for i in range(0,6):
            selected_chrs = random.choices(host_chrs, k=1)
            if i > 2:
                selected_chrs = random.choices(host_chrs, k=2)
                while len(set(selected_chrs)) != 2:
                    selected_chrs = random.choices(host_chrs, k=2)
            if i >= 5:
                selected_chrs = random.choices(host_chrs, k=3)
                while len(set(selected_chrs)) != 3:
                    selected_chrs = random.choices(host_chrs, k=3)
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
            run_local(out_dir,script_root,vc,v_len,selected_chrs)
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
                i_start = random.randint(500,v_len - 1500)
                i_end = random.randint(i_start + 500,v_len)
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