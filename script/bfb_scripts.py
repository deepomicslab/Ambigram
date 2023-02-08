import argparse
from cmath import inf
import sys
import os

class MainArgParser:
    def __init__(self):
        parser = argparse.ArgumentParser(prog='preBFB')
        parser.add_argument(dest='subfunc', help='Available sub-functions: cluster_sv, generate_seg, call_depth, generate_lh')
        args = parser.parse_args(sys.argv[1:2])
        getattr(self, args.subfunc)()

    """ Simulate reads with different simulators and protocols """
    def seg2bed(self):
        parser = argparse.ArgumentParser(description='Given seg file, generate .bed file.')
        parser.add_argument('-s', '--seg_file', dest='seg', required=True, help='Path to seg file')
        parser.add_argument('-p', '--prefix', dest='prefix', required=False,  default = 'test', help='Prefix of output .bed file')
        
        args = parser.parse_args(sys.argv[2:])
        segments = []
        for line in open(args.seg, 'r').readlines():
            info = line.strip('\n').split(' ')
            chrName, bkp = info[0].split(':')
            start, end = bkp.split('-')
            segments.append([chrName, start, end]+info[1:])
        with open('{}_seg.bed'.format(args.prefix), 'w') as output:
            for segment in segments:
                output.write('\t'.join(segment)+'\n')
            output.close()

    def getFasta(self):
        parser = argparse.ArgumentParser(description='Given .bed file, generate .fa file.')
        parser.add_argument('-b', '--bed_file', dest='bed', required=True, help='Path to .bed file')
        parser.add_argument('-r', '--reference', dest='ref', required=True, help='Path to reference genome file (.fa)')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False,  default = 'test', help='Sample name for output file')
        
        args = parser.parse_args(sys.argv[2:])
        # generate .fa
        cmd = 'bedtools getfasta -fi {} -bed {} -s -fo {}.fa'.format(args.ref, args.bed, args.sampleName)
        os.system(cmd)
        # merge short sequences
        res = ">BFB\n"
        for line in open('{}.fa'.format(args.sampleName), "r").readlines():
            if line[0] == '>':
                continue
            res += line.strip('\n')
        fasta = open('{}.fa'.format(args.sampleName), "w")
        fasta.write(res)
        fasta.close()

    def simulate_PE(slef):
        parser = argparse.ArgumentParser(description='Simulate pair-end reads.')
        parser.add_argument('-f', '--fasta_file', dest='fasta', required=True, help='Path to the input .fa file')
        parser.add_argument('-r', '--reference', dest='ref', required=True, help='Path to reference genome file (.fa)')
        parser.add_argument('-p', '--tumor_purity', dest='purity', required=False, type=float, default = 1, help='Sample tumor purity (default: 1)')
        parser.add_argument('-n', '--normal_bam', dest='normalBam', required=False, default=None, help='Path to the normal .bam file')
        parser.add_argument('-c', '--coverage', dest='coverage', required=False, type=int, default=30, help='Average coverage on each base (default: 30)')
        parser.add_argument('-l', '--read_length', dest='length', required=False, type=int, default=150, help='Length of each read (default: 150)')
        parser.add_argument('-i', '--insertion_size', dest='insertion', required=False, type=int, default=350, help='Insertion distance between the paried reads (default: 350)')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False,  default = 'test', help='Sample name for output file')
        
        args = parser.parse_args(sys.argv[2:])
        
        # simulate reads
        fa_length = 0
        for line in open(args.fasta, 'r').readlines():
            if line.startswith('>'):
                continue
            fa_length += len(line.strip('\n'))
        read_length = fa_length*args.coverage/args.length/2
        cmd = 'wgsim -1 {} -2 {} -d {} -N {} -r 0 -R 0 -X 0 -e 0 {} test_r1.fq test_r2.fq'\
            .format(args.length, args.length, args.insertion, read_length, args.fasta)
        os.system(cmd)
        # align reads with the reference
        cmd = 'bwa mem -t 8 -o {}.sam {} test_r1.fq test_r2.fq && samtools sort -O BAM {}.sam -o {}.bam --threads 8 && samtools index {}.bam'\
            .format(args.sampleName, args.ref, args.sampleName, args.sampleName, args.sampleName)
        os.system(cmd)
        if args.purity < 1:
            cmd = 'samtools view -b -s {} {} > temp1.bam && samtools view -b -s {} {}.bam > temp2.bam'\
                .format(1-args.purity, args.normalBam, args.purity, args.sampleName)
            cmd += ' && samtools merge {}.bam temp1.bam temp2.bam -f && samtools index {}.bam'\
                .format(args.sampleName, args.sampleName)
            cmd += ' && rm -f temp1.bam && rm -f temp2.bam'
            print(cmd)
            os.system(cmd)
        # call SV
        cmd = 'svaba run -t {}.bam -G {} -a {}'.format(args.sampleName, args.ref, 'sample')
        os.system(cmd)
        # remove .sam and .bai
        cmd = 'rm -f {}.sam {}.bam.bai'.format(args.sampleName, args.sampleName)
        os.system(cmd)

    def simulate_PB(slef):
        parser = argparse.ArgumentParser(description='Simulate PacBio long reads.')
        parser.add_argument('-f', '--fasta_file', dest='fasta', required=True, help='Path to the input .fa file')
        parser.add_argument('-r', '--reference', dest='ref', required=True, help='Path to reference genome file (.fa)')
        parser.add_argument('-p', '--tumor_purity', dest='purity', required=False, type=float, default = 1, help='Sample tumor purity (default: 1)')
        parser.add_argument('-n', '--normal_bam', dest='normalBam', required=False, default=None, help='Path to the normal .bam file')
        parser.add_argument('-c', '--coverage', dest='coverage', required=False, type=int, default=30, help='Average coverage on each base (default: 30)')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False,  default = 'test', help='Sample name for output file')
        
        args = parser.parse_args(sys.argv[2:])
        
        # simulate reads
        cmd = '~/pbsim3/build/src/pbsim --strategy wgs --method qshmm --qshmm ~/pbsim3/data/QSHMM-RSII.model --depth {} --genome {} --prefix sample_pb'\
            .format(args.coverage, args.fasta)
        # cmd = '~/pbsim2/build/src/pbsim --prefix sample_pb --depth {} --hmm_model ~/pbsim2/data/P6C4.model {}'\
        #     .format(args.coverage, args.fasta)
        os.system(cmd)
        # align reads with the reference
        cmd = 'ngmlr -t 8 -r ~/csv/reference/hg38.fa -q sample_pb_0001.fastq -o {}.sam && samtools sort -O BAM {}.sam -o {}.bam --threads 8 && samtools index {}.bam'\
            .format(args.sampleName, args.sampleName, args.sampleName, args.sampleName)
        os.system(cmd)
        if args.purity < 1:
            cmd = 'samtools view -b -s {} {} > temp1.bam && samtools view -b -s {} {}.bam > temp2.bam'\
                .format(1-args.purity, args.normalBam, args.purity, args.sampleName)
            cmd += ' && samtools merge {}.bam temp1.bam temp2.bam -f && samtools index {}.bam'\
                .format(args.sampleName, args.sampleName)
            cmd += ' && rm -f temp1.bam && rm -f temp2.bam'
            print(cmd)
            os.system(cmd)
        # call SV
        cmd = 'sniffles -m {}.bam -v sample_pb.vcf -s 1'.format(args.sampleName)
        os.system(cmd)
        # remove .sam and .bai
        cmd = 'rm -f {}.sam {}.bam.bai'.format(args.sampleName, args.sampleName)
        os.system(cmd)

    def simulate_ONT(slef):
        parser = argparse.ArgumentParser(description='Simulate ONT long reads.')
        parser.add_argument('-f', '--fasta_file', dest='fasta', required=True, help='Path to the input .fa file')
        parser.add_argument('-r', '--reference', dest='ref', required=True, help='Path to reference genome file (.fa)')
        parser.add_argument('-p', '--tumor_purity', dest='purity', required=False, type=float, default = 1, help='Sample tumor purity (default: 1)')
        parser.add_argument('-n', '--normal_bam', dest='normalBam', required=False, default=None, help='Path to the normal .bam file')
        parser.add_argument('-c', '--coverage', dest='coverage', required=False, type=int, default=30, help='Average coverage on each base (default: 30)')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False,  default = 'test', help='Sample name for output file')
        
        args = parser.parse_args(sys.argv[2:])
        
        # simulate reads
        cmd = '~/pbsim3/build/src/pbsim --strategy wgs --method qshmm --qshmm ~/pbsim3/data/QSHMM-ONT.model --depth {} --genome {} --prefix sample_ont'\
            .format(args.coverage, args.fasta)
        # cmd = '~/pbsim2/build/src/pbsim --prefix sample_ont --depth {} --hmm_model ~/pbsim2/data/R103.model {}'\
        #     .format(args.coverage)
        os.system(cmd)
        # align reads with the reference
        cmd = 'ngmlr -t 8 -r ~/csv/reference/hg38.fa -q sample_ont_0001.fastq -o {}.sam && samtools sort -O BAM {}.sam -o {}.bam --threads 8 && samtools index {}.bam'\
            .format(args.sampleName, args.sampleName, args.sampleName, args.sampleName)
        os.system(cmd)
        if args.purity < 1:
            print('Merge the normal and tumor samples ...')
            cmd = 'samtools view -b -s {} {} > temp1.bam && samtools view -b -s {} {}.bam > temp2.bam'\
                .format(1-args.purity, args.normalBam, args.purity, args.sampleName)
            cmd += ' && samtools merge {}.bam temp1.bam temp2.bam -f && samtools index {}.bam'\
                .format(args.sampleName, args.sampleName)
            cmd += ' && rm -f temp1.bam && rm -f temp2.bam'
            print(cmd)
            os.system(cmd)
        # call SV
        cmd = 'sniffles -m {}.bam -v sample_ont.vcf -s 1'.format(args.sampleName)
        os.system(cmd)
        # remove .sam and .bai
        cmd = 'rm -f {}.sam {}.bam.bai'.format(args.sampleName, args.sampleName)
        os.system(cmd)

    def simulate_10x(slef):
        parser = argparse.ArgumentParser(description='Simulate 10x linked reads.')
        parser.add_argument('-f', '--fasta_file', dest='fasta', required=True, help='Path to the input .fa file')
        parser.add_argument('-r', '--reference', dest='ref', required=True, help='Path to reference genome file (.fa)')
        parser.add_argument('-p', '--tumor_purity', dest='purity', required=False, type=float, default = 1, help='Sample tumor purity (default: 1)')
        parser.add_argument('-t', '--tumor_bam', dest='tumorBam', required=False, default=None, help='Path to the tumor .bam file')
        parser.add_argument('-n', '--normal_bam', dest='normalBam', required=False, default=None, help='Path to the normal .bam file')
        parser.add_argument('-c', '--coverage', dest='coverage', required=False, type=int, default=30, help='Average coverage on each base (default: 30)')
        parser.add_argument('-l', '--read_length', dest='length', required=False, type=int, default=150, help='Length of each read (default: 150)')
        parser.add_argument('-i', '--insertion_size', dest='insertion', required=False, type=int, default=350, help='Insertion distance between the paried reads (default: 350)')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False,  default = 'test', help='Sample name for output file')
        
        args = parser.parse_args(sys.argv[2:])
        
        # simulate reads
        fa_length = 0
        for line in open(args.fasta, 'r').readlines():
            if line.startswith('>'):
                continue
            fa_length += len(line.strip('\n'))
        read_pairs = fa_length*args.coverage/(args.length*2)/1000000
        pe_length =  fa_length*args.coverage/2
        partition_lower = (pe_length/50/1000+1)/1000
        cmd = '~/LRSIM/simulateLinkedReads.pl -g {} -p ./10x_fastq/{} -d 1 -1 {} -i {} -x {} -f 50 -t {} -m 1 -o -4 1 -7 1 -e 0'\
            .format(args.fasta, args.sampleName, fa_length, args.insertion, read_pairs, partition_lower)
        os.system(cmd)
        
        # alignment
        cmd = 'bash ./10x_fastq/simulate_10x.sh tenx_simulation 30 150 350 {} {} {}'\
            .format(args.sampleName, args.fasta, args.ref)
        os.system(cmd)
        # call SV
        if args.purity < 1:
            print('Merge the normal and tumor samples ...')
            cmd = 'samtools view -b -s {} {} > temp1.bam && samtools view -b -s {} {} > temp2.bam'\
                .format(1-args.purity, args.normalBam, args.purity, args.tumorBam)
            cmd += ' && samtools merge {}.bam temp1.bam temp2.bam -f && samtools index {}.bam'\
                .format(args.sampleName, args.sampleName)
            cmd += ' && rm -f temp1.bam && rm -f temp2.bam'
            print(cmd)
            os.system(cmd)
        cmd = 'svaba run -t {}.bam -G {} -a {}'.format(args.sampleName, args.ref, 'sample')
        os.system(cmd)

    def sniffles2sv(self):
        parser = argparse.ArgumentParser(description='Extract SV information from sniffles VCF file')
        parser.add_argument('-v', '--vcf', dest='vcf', required=True, help='Input VCF file')
        parser.add_argument('-p', '--prefix', dest='prefix', required=False, default = 'test', help='Prefix of output file')
        args = parser.parse_args(sys.argv[2:])

        sv = []
        for line in open(args.vcf, "r").readlines():
            if line.startswith('#'):
                continue
            info = line.strip('\n').split('\t')
            prop = {}
            for elem in info[7].split(';')[1:]:
                prop[elem.split('=')[0]] = elem.split('=')[1]
            chr1, pos1 = info[0], info[1]
            chr2, pos2 = prop['CHR2'], prop['END']
            str1, str2 = '+', '+'
            if prop['STRANDS'] == '++':
                str1, str2 = '+', '-'
            elif prop['STRANDS'] == '--':
                str1, str2 = '-', '+'
            elif prop['STRANDS'] == '+-':
                str1, str2 = '+', '+'
            elif prop['STRANDS'] == '-+':
                str1, str2 = '-', '-'
            depth = info[-1].split(":")[-1]
            sv.append([chr1, pos1, str1, chr2, pos2, str2, depth])
        print(sv)
        with open('{}_sv.txt'.format(args.prefix), 'w') as sv_file:
            sv_file.write('chr_3p	bkp_3p	str_3p	chr_5p	bkp_5p	str_5p	depth\n')
            for info in sv:
                sv_file.write('\t'.join(info)+'\n')
            sv_file.close()

    def svaba2sv(self):
        parser = argparse.ArgumentParser(description='Extract SV information from svaba VCF file')
        parser.add_argument('-v', '--vcf', dest='vcf', required=True, help='Input VCF file')
        parser.add_argument('-p', '--prefix', dest='prefix', required=False, default = 'test', help='Prefix of output file')
        args = parser.parse_args(sys.argv[2:])

        sv = []
        for line in open(args.vcf, "r").readlines():
            if line.startswith('#'):
                continue
            info = line.strip('\n').split('\t')
            if info[2][-1] == '2':
                continue
            prop = info[7].split(';')
            end = info[4].split('[')
            str1, str2 = '+', '+'
            if ']' in info[4]:
                end = info[4].split(']')
                str2 = '-'
            if end[0] == '':
                str1 = '-'
            chr1, bkp1 = info[0], info[1]
            chr2, bkp2 = end[1].split(':')
            key, num = info[8].split(':'), info[12].split(':')
            data = {}
            for i in range(len(key)):
                data[key[i]] = num[i]
            ad, dp = int(data['AD']), int(data['DP'])
            sv.append([chr1, bkp1, str1, chr2, bkp2, str2, data['AD']])
            id = '{}:{}:{}:{}:{}:{}'.format(chr1, bkp1, str1, chr2, bkp2, str2)
        with open('{}_sv.txt'.format(args.prefix), 'w') as sv_file:
            sv_file.write('chr_3p	bkp_3p	str_3p	chr_5p	bkp_5p	str_5p	depth\n')
            for info in sv:
                sv_file.write('\t'.join(info)+'\n')
            sv_file.close()

    def OM2juncs(self):
        parser = argparse.ArgumentParser(description='Convert optical mapping alignment from SegAligner to .junc file')
        parser.add_argument('-i', '--input', dest='input', required=True, help='Input OM alignment file from SegAligner')
        parser.add_argument('-p', '--prefix', dest='prefix', required=False, default = 'test', help='Prefix of output file')
        args = parser.parse_args(sys.argv[2:])

        res = ''
        for line in open(args.input, 'r').readlines():
            if line.startswith('#'):
                continue
            seg = line.split('\t')[0]
            print(seg)
            if seg.startswith('-'):
                res += seg[1:]+seg[0]+' '
            else:
                res += seg+'+ '
        with open('{}.juncs'.format(args.prefix), 'w') as output:
            output.write(res[:-1])
            output.close()

    """ Cluster complex SV information """
    def compare(self, sv1, sv2):
        if int(sv1[1]) < int(sv2[1]):
            return -1
        elif int(sv1[1]) > int(sv2[1]):
            return 1
        if int(sv1[4]) < int(sv2[4]):
            return -1
        elif int(sv1[4]) > int(sv2[4]):
            return 1
        return 0

    def min_dis(self, sv1, sv2):
        diff1 = abs(int(sv1[1])-int(sv2[1]))
        if sv1[0] != sv2[0]:
            diff1 = inf
        diff2 = abs(int(sv1[1])-int(sv2[4]))
        if sv1[0] != sv2[3]:
            diff2 = inf
        diff3 = abs(int(sv1[4])-int(sv2[1]))
        if sv1[3] != sv2[0]:
            diff3 = inf
        diff4 = abs(int(sv1[4])-int(sv2[4]))
        if sv1[3] != sv2[3]:
            diff4 = inf
        return min(diff1,diff2,diff3,diff4)
    
    def setRange(self, chr_range, sv):
        if sv[0] in chr_range.keys():
            chr_range[sv[0]][0] = min(chr_range[sv[0]][0], int(sv[1]))
            chr_range[sv[0]][1] = max(chr_range[sv[0]][1], int(sv[1]))
        else:
            chr_range[sv[0]] = [int(sv[1]), int(sv[1])]
        if sv[3] in chr_range.keys():
            chr_range[sv[3]][0] = min(chr_range[sv[3]][0], int(sv[4]))
            chr_range[sv[3]][1] = max(chr_range[sv[3]][1], int(sv[4]))
        else:
            chr_range[sv[3]] = [int(sv[4]), int(sv[4])]
        return chr_range
    
    def check_range(self, chr_range, max_range):
        for val in chr_range.values():
            if val[1]-val[0] > max_range:
                return False
        return True

    def hasFBI(self, sv_id, sv):
        for i in sv_id:
            if sv[i][0] == sv[i][3] and sv[i][2] != sv[i][5]:
                return True
        return False

    def cluster_sv(self):
        parser = argparse.ArgumentParser(description='Cluster complex SV based on breakpoint distance.')
        parser.add_argument('-sv', '--sv_file', dest='svPath', required=True, help='Path to SV file')
        parser.add_argument('-d', '--max_dis', dest='maxDis', required=False, type=int, default=1000000, help='Maximum distance of two SVs grouped in a cluster')
        parser.add_argument('-r', '--max_range', dest='maxRange', required=False, type=int, default=10000000, help='Maximum range of a cluster')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False,  default = 'test', help='Sample name for output file')
        args = parser.parse_args(sys.argv[2:])
        # read all SVs
        juncs = []
        for line in open(args.svPath).readlines()[1:]:
            if line == '\n':
                continue
            info = line.strip('\n').split('\t')
            if info[2] == '-' and info[5] == '-':
                info[2], info[5] = '+', '+'
                info[0], info[3] = info[3], info[0]
                info[1], info[4] = info[4], info[1]
            juncs.append(info)
        # sort juncs in chromosome
        juncs = sorted(juncs, key=lambda x: (x[0], self.compare))
        # cluster SV junctions based on distance
        cluster = []
        svIdx = list(range(0, len(juncs))) # sv index for selection
        while len(svIdx) > 0:
            subcluster, queue = [svIdx[0]], [svIdx[0]] # index of sv_info
            sv, chr_range = juncs[svIdx[0]], {}
            self.setRange(chr_range, sv)
            svIdx.pop(0)
            while len(queue) > 0:
                idx = queue[0]
                queue.pop(0)
                for i in svIdx:
                    if self.min_dis(juncs[i], juncs[idx]) < args.maxDis:
                        temp_range = chr_range.copy()
                        self.setRange(temp_range, juncs[i])
                        if self.check_range(temp_range, args.maxRange):
                            self.setRange(chr_range, juncs[i])
                            queue.append(i)
                            subcluster.append(i)
                            svIdx.remove(i)
                            print(subcluster)
            if self.hasFBI(subcluster, juncs) == True:
                cluster.append(subcluster)
        # output a set of sv.txt files
        for i in range(len(cluster)):
            svFile = open('{}_{}_sv.txt'.format(args.sampleName, i+1), "w")
            res = 'chrom_5p\tbkpos_5p\tstrand_5p\tchrom_3p\tbkpos_3p\tstrand_3p\tavg_cn\n'
            for idx in cluster[i]:                
                res += '\t'.join(juncs[idx])+'\n'
            svFile.write(res)
            svFile.close()

    """ Generate SEG file from SV and BAM files """
    # convert coverage depth to copy number
    def depth2cn(self, sampleDepth, wgsDepth, purity):
        ploidy = 2
        haploDepth = (wgsDepth*purity/ploidy)
        return sampleDepth/haploDepth

    # generate segments based on SV breakpoints
    def generate_seg(self):
        import pysam
        parser = argparse.ArgumentParser(description='Generate SEG file from SV and BAM files.')
        parser.add_argument('-sv', '--sv_file', dest='svPath', required=True, help='Path to SV file')
        parser.add_argument('-bam', '--bam_file', dest='bamPath', required=False, help='Path to BAM file')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False,  default = 'test', help='Sample name for output file')
        parser.add_argument('-d', '--wgs_depth', dest='wgsDepth', required=False, type=int, default = 30, help='The whole genome average depth (default: 30)')
        parser.add_argument('-p', '--tumor_purity', dest='purity', required=False, type=float, default = 1, help='Sample tumor purity (default: 1)')
        args = parser.parse_args(sys.argv[2:])
        # find all breakpoints on each chromosome
        sv, pos = [], {}
        for line in open(args.svPath, 'r').readlines()[1:]:
            info = line.strip('\n').split('\t')
            info[1], info[4] = int(info[1]), int(info[4])
            sv.append(info)
            if info[0] in pos.keys():
                pos[info[0]].append(info[1])
            else:
                pos[info[0]] = [info[1]]
            if info[3] in pos.keys():
                pos[info[3]].append(info[4])
            else:
                pos[info[3]] = [info[4]]
        for key, arr in pos.items():
            pos[key] = list(set(arr))
            pos[key].sort()
            pos[key].insert(0, max(1, pos[key][0]-1000))
            pos[key].append(pos[key][-1]+1000)

        # call coverage depth of segments from BAM
        segDepth = {}
        if args.bamPath == None:
            for chr, arr in pos.items():
                for n in range(1, len(arr)):
                    key = chr+':'+str(arr[n-1])+'-'+str(arr[n])
                    segDepth[key] = 100
        else:
            bam = pysam.AlignmentFile(args.bamPath, 'rb')
            for key, arr in pos.items():
                for n in range(1, len(arr)):
                    cnt = bam.count_coverage(key, arr[n-1], arr[n], quality_threshold = 0)
                    posDepth = []
                    for i in range(len(cnt[0])):
                        temp = 0
                        for j in range(4):
                            temp += cnt[j][i]
                        posDepth.append(temp)
                    name = key+':'+str(arr[n-1])+'-'+str(arr[n])
                    segDepth[name] = sum(posDepth)/len(posDepth) # average on read depth of all positions
        # output seg.txt file
        if args.wgsDepth!=30 and args.purity!=1:
            for key, value in segDepth.items():
                segDepth[key] = self.depth2cn(value, args.wgsDepth, args.purity)
                print(key, segDepth[key])
        resStr = ''
        for key, value in segDepth.items():
            resStr += key+'\t'+str(value)+'\n'
        segFile = open('{}_seg.txt'.format(args.sampleName), 'w')
        segFile.write(resStr)
        segFile.close()

    """ call base depth from BAM """
    def call_depth(self):
        import pysam
        parser = argparse.ArgumentParser(description='Extract coverage depth of bases on regions.')
        parser.add_argument('-seg', '--seg_file', dest='segPath', required=True, help='Path to SEG file')
        parser.add_argument('-bam', '--bam_file', dest='bamPath', required=True, help='Path to BAM file')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False,  default = 'test', help='Sample name for output file')
        args = parser.parse_args(sys.argv[2:])
        # call coverage depth
        bam = pysam.AlignmentFile(args.bamPath, 'rb')
        res = ''
        for line in open(args.segPath).readlines():
            line = line.strip('\n').split('\t')[0]
            ref, bkp = line.split(':')[0], line.split(':')[1].split('-')        
            cnt = bam.count_coverage(ref, int(bkp[0]), int(bkp[1])+1, quality_threshold = 0)
            pos = int(bkp[0])
            for i in range(len(cnt[0])):
                temp = 0
                for j in range(4):
                    temp += cnt[j][i]
                res += ref+'\t'+str(pos)+'\t'+str(temp)+'\n'
                pos = pos + 1
        coverageFile = open('{}_coverage.txt'.format(args.sampleName), 'w')
        coverageFile.write(res)
        coverageFile.close()
    
    """ Generate .lh file by integrating sv.txt and seg.txt """
    def findSegment(self, segs, bkp, isStart): # bkp = [chr, pos, strand]
        isLeft = True # left breakpoint
        if (isStart == True and bkp[2] == '+') or (isStart == False and bkp[2] == '-'):
            isLeft = False # right breakpoint
        segID = len(segs)
        minDis = float('inf')
        for seg in segs:
            if bkp[0] == seg[1]:
                if isLeft==True and abs(int(seg[2])-int(bkp[1])) < minDis:
                    segID, minDis = seg[0], abs(int(seg[2])-int(bkp[1]))
                elif isLeft==False and abs(int(seg[3])-int(bkp[1])) < minDis:
                    segID, minDis = seg[0], abs(int(seg[3])-int(bkp[1]))
        return segID
    
    def hasDuplicateSV(self, sv, info): # info = [seg1, str1, seg2, str2]
        for junc in sv:
            if junc[0] == info[0] and junc[2] == info[2]:
                if junc[1] == info[1] and junc[3] == info[3]:
                    return sv.index(junc)
            elif junc[0] == info[2] and junc[2] == info[0]:
                if info[1] != info[3]:
                    if junc[1] == info[1] and junc[3] == info[3]:
                        return sv.index(junc)
                else:
                    if junc[1] != info[1] and junc[3] != info[3]:
                        return sv.index(junc)
        return -1

    def generate_lh(self):
        parser = argparse.ArgumentParser(description='Generate .lh file by integrating sv.txt and seg.txt.')
        parser.add_argument('-sv', '--sv_file', dest='svPath', required=True, help='Path to SV file')
        parser.add_argument('-seg', '--seg_file', dest='segPath', required=True, help='Path to SEG file')
        parser.add_argument('-c', '--coverage', dest='coverage', required=False, type=int, default=30, help='Average coverage on each base (default: 30)')
        parser.add_argument('-p', '--tumor_purity', dest='purity', required=False, type=float, default = 1, help='Sample tumor purity (default: 1)')
        parser.add_argument('-d', '--is_depth', dest='isDepth', required=False, default = False, help='Indicate either coverage depth or copy number is provided')
        parser.add_argument('-d1', '--is_seg_depth', dest='isSegDepth', required=False, default = False, help='Indicate either segment coverage depth or copy number is provided')
        parser.add_argument('-d2', '--is_sv_depth', dest='isSVDepth', required=False, default = False, help='Indicate either SV coverage depth or copy number is provided')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False,  default = 'test', help='Sample name for output file')
        parser.add_argument('-pr', '--property', dest='prop', required=False,  default = '', help='PROP inforamtion in .lh file')

        args = parser.parse_args(sys.argv[2:])
        # read segments and construct normal junctions
        segs = []
        sourceSegs, sinkSegs = [1], []
        cnt = 1
        for line in open(args.segPath, 'r').readlines():
            info = line.strip('\n').split('\t')
            [chrName, interval] = info[0].split(':')
            segs.append([cnt, chrName, interval.split('-')[0], interval.split('-')[1], info[1]])
            # print('seg{} length:{} mid:{}'.format(cnt, int(segs[-1][3])-int(segs[-1][2]), (int(segs[-1][3])+int(segs[-1][2]))/2))
            if chrName != segs[sourceSegs[-1]-1][1]:
                sinkSegs.append(cnt-1)
                sourceSegs.append(cnt)
            cnt += 1
        sinkSegs.append(cnt-1)
        # read SVs
        sv = []
        for line in open(args.svPath, 'r').readlines()[1:]:
            info = line.strip('\n').split('\t')
            # if info[0] == info[3] and info[2] == info[5]:
            #     continue
            seg1, seg2 = self.findSegment(segs, info[:3], True), self.findSegment(segs, info[3:6], False)
            if int(seg1)+1 == int(seg2) and info[2] == info[5] and info[2] == '+':
                continue
            if int(seg1) == int(seg2)+1 and info[2] == info[5] and info[2] == '-':
                continue
            # if info[0] == info[3] and abs(seg1-seg2) > 2:
            #     continue
            junc_index = self.hasDuplicateSV(sv, [seg1, info[2], seg2, info[5]])
            if junc_index != -1:
                if float(info[6]) > float(sv[junc_index][-1]):
                    sv[junc_index][-1] = info[6]
            else:
                sv.append([seg1, info[2], seg2, info[5], info[6]])
            
        # output .lh file
        res = '''SAMPLE group1
AVG_CHR_SEG_DP {}
AVG_WHOLE_HOST_DP {}
AVG_JUNC_DP {}
PURITY {}
AVG_TUMOR_PLOIDY 2
PLOIDY 2m1
VIRUS_START 7
SOURCE {}
SINK {}
'''.format(args.coverage, args.coverage, args.coverage, args.purity, ','.join(str(e) for e in sourceSegs), ','.join(str(e) for e in sinkSegs))
        for i in range(len(segs)):
            if segs[i][1] == 'chr18':
                segs[i][1] = 'virus'
        if args.isSegDepth == False and args.isDepth == False:# segment
            for seg in segs:
                res += 'SEG H:{}:{}:{}:{} {} {}\n'.format(seg[0], seg[1], seg[2], seg[3], float(seg[4])*30, seg[4])
        else:
            for seg in segs:
                res += 'SEG H:{}:{}:{}:{} {} {}\n'.format(seg[0], seg[1], seg[2], seg[3], seg[4], -1)
        if args.isSVDepth == False and args.isDepth == False:# SV junction
            for junc in sv:
                res += 'JUNC H:{}:{} H:{}:{} {} {} U B\n'.format(junc[0], junc[1], junc[2], junc[3], float(junc[4])*30, junc[4])
        else:
            for junc in sv:
                res += 'JUNC H:{}:{} H:{}:{} {} {} U B\n'.format(junc[0], junc[1], junc[2], junc[3], junc[4], -1)
        res += args.prop
        # if args.isSegDepth == False and args.isDepth == False:# normal junction
        #     for junc in normal:
        #         res += 'JUNC H:{}:{} H:{}:{} {} {} U B\n'.format(junc[0], junc[1], junc[2], junc[3], float(junc[4])*15, junc[4])
        # else:
        #     for junc in normal:
        #         res += 'JUNC H:{}:{} H:{}:{} {} {} U B\n'.format(junc[0], junc[1], junc[2], junc[3], junc[4], -1)
        lhFile = open('{}.lh'.format(args.sampleName), 'w')
        lhFile.write(res)
        lhFile.close()

if __name__ == "__main__":
    MainArgParser()