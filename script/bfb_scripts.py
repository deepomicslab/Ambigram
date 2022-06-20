import argparse
from cmath import inf
import sys

class MainArgParser:
    def __init__(self):
        parser = argparse.ArgumentParser(prog='preBFB')
        parser.add_argument(dest='subfunc', help='Available sub-functions: cluster_sv, generate_seg, call_depth, generate_lh')
        args = parser.parse_args(sys.argv[1:2])
        getattr(self, args.subfunc)()

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

    def cluster_sv(self):
        import functools
        parser = argparse.ArgumentParser(description='Cluster complex SV based on breakpoint distance.')
        parser.add_argument('-sv', '--sv_file', dest='svPath', required=True, help='Path to SV file')
        parser.add_argument('-d', '--max_dis', dest='maxDis', required=False, default=1000000, help='Maximum distance of two SVs in a cluster')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False, default = 'sample', help='Sample name for output file')
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
            subcluster = [svIdx[0]] # index of sv_info
            queue = [svIdx[0]]
            svIdx.pop(0)
            while len(queue) > 0:
                idx = queue[0]
                queue.pop(0)
                for i in svIdx:
                    if self.min_dis(juncs[idx], juncs[i]) > int(args.maxDis):
                        continue
                    queue.append(i)
                    subcluster.append(i)
                    svIdx.remove(i)
            cluster.append(subcluster)
        # output a set of sv.txt files
        for i in range(len(cluster)):
            svFile = open('{}{}_sv.txt'.format(args.sampleName, i+1), "w")
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
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False, default = 'sample', help='Sample name for output file')
        parser.add_argument('-d', '--wgs_depth', dest='wgsDepth', required=False, type=int, default = 100, help='The whole genome average depth (default: 100)')
        parser.add_argument('-p', '--tumor_purity', dest='purity', required=False, type=int, default = 1, help='Sample tumor purity (default: 1)')
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
            pos[key].insert(0, pos[key][0]-1000)
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
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False, default = 'sample', help='Sample name for output file')
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
        idx = 2 # left breakpoint
        if (isStart == True and bkp[2] == '+') or (isStart == False and bkp[2] == '-'):
            idx = 3 # right breakpoint
        segID = 1
        for seg in segs:
            if bkp[0] == seg[1]:
                if idx == 2 and int(bkp[1]) < int(seg[idx]):
                    segID = seg[0]-1
                    break
                if idx == 3 and int(bkp[1]) <= int(seg[idx]):
                    segID = seg[0]
                    break
        return segID

    def generate_lh(self):
        parser = argparse.ArgumentParser(description='Generate .lh file by integrating sv.txt and seg.txt.')
        parser.add_argument('-sv', '--sv_file', dest='svPath', required=True, help='Path to SV file')
        parser.add_argument('-seg', '--seg_file', dest='segPath', required=True, help='Path to SEG file')
        parser.add_argument('-s', '--sample_name', dest='sampleName', required=False, default = 'sample', help='Sample name for output file')

        args = parser.parse_args(sys.argv[2:])
        # read segments and construct normal junctions
        segs, normal = [], []
        sourceSegs, sinkSegs = [1], []
        cnt = 1
        for line in open(args.segPath, 'r').readlines():
            info = line.strip('\n').split('\t')
            [chrName, interval] = info[0].split(':')
            segs.append([cnt, chrName, interval.split('-')[0], interval.split('-')[1], info[1]])
            if chrName != segs[sourceSegs[-1]-1][1]:
                sinkSegs.append(cnt-1)
                sourceSegs.append(cnt)
            elif cnt > 1:
                cn = (float(segs[cnt-2][-1])+float(segs[cnt-1][-1]))/2
                normal.append([cnt-1, '+', cnt, '+', str(cn)])
            cnt += 1
        sinkSegs.append(cnt-1)
        # read SVs
        sv = []
        for line in open(args.svPath, 'r').readlines()[1:]:
            info = line.strip('\n').split('\t')
            if info[0] == info[3] and info[2] == info[5]:
                continue
            seg1, seg2 = self.findSegment(segs, info[:3], True), self.findSegment(segs, info[3:6], False)
            if info[0] == info[3] and abs(seg1-seg2) > 2:
                continue
            sv.append([seg1, info[2], seg2, info[5], info[6]])
            
        # output .lh file
        res = '''SAMPLE group1
AVG_CHR_SEG_DP 30
AVG_WHOLE_HOST_DP 30
AVG_JUNC_DP 30
PURITY 1
AVG_TUMOR_PLOIDY 2
PLOIDY 2m1
VIRUS_START 7
SOURCE {}
SINK {}
'''.format(','.join(str(e) for e in sourceSegs), ','.join(str(e) for e in sinkSegs))   
        for seg in segs:
            res += 'SEG H:{}:{}:{}:{} {} {}\n'.format(seg[0], seg[1], seg[2], seg[3], float(seg[4])*15, seg[4])
        for junc in sv:
            res += 'JUNC H:{}:{} H:{}:{} {} {} U B\n'.format(junc[0], junc[1], junc[2], junc[3], float(junc[4])*15, junc[4])
        for junc in normal:
            res += 'JUNC H:{}:{} H:{}:{} {} {} U B\n'.format(junc[0], junc[1], junc[2], junc[3], float(junc[4])*15, junc[4])
        lhFile = open('{}.lh'.format(args.sampleName), 'w')
        lhFile.write(res)
        lhFile.close()

if __name__ == "__main__":
    MainArgParser()