import pysam
import argparse

def getSV(svPath):
    svFile = open(svPath, 'r')
    sv, pos = [], {}
    for line in svFile.readlines()[1:]:
        info = line.strip('\n').split('\t')
        info[1] = int(info[1])
        info[4] = int(info[4])
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
    # print(sv, pos)
    return sv, pos

def getSegDepth(bamPath, pos):
    bam = pysam.AlignmentFile(bamPath, 'rb')
    depth = {}
    for key, arr in pos.items():
        for n in range(1, len(arr)):
            cnt = bam.count_coverage(key, arr[n-1]-1, arr[n], quality_threshold = 0)
            posDepth = []
            for i in range(len(cnt[0])):
                temp = 0
                for j in range(4):
                    temp += cnt[j][i]
                posDepth.append(temp)
            name = key+':'+str(arr[n-1])+'-'+str(arr[n])
            depth[name] = sum(posDepth)/len(posDepth) # average on read depth of all positions        
    return depth

def depth2cn(sampleDepth, wgsDepth, purity):
    ploidy = 2
    haploDepth = (wgsDepth*purity/ploidy)
    return sampleDepth/haploDepth


def main():
    parser = argparse.ArgumentParser(description='Generate SEG file from SV and BAM files.')
    parser.add_argument('-sv', '--sv_file', dest='svPath', required=True, help='Path to SV file')
    parser.add_argument('-bam', '--bam_file', dest='bamPath', required=True, help='Path to BAM file')
    parser.add_argument('-s', '--sample_name', dest='sampleName', required=False, default = 'sample', help='Sample name for output file')
    parser.add_argument('-d', '--wgs_depth', dest='wgsDepth', required=False, type=int, default = 100, help='The whole genome average depth (default: 100)')
    parser.add_argument('-p', '--tumor_purity', dest='purity', required=False, type=int, default = 1, help='Sample tumor purity (default: 1)')
    args = parser.parse_args()
    sv, pos = getSV(args.svPath)
    segDepth = getSegDepth(args.bamPath, pos)
    for key, value in segDepth.items():
        segDepth[key] = depth2cn(value, args.wgsDepth, args.purity)
        print(key, segDepth[key])
    resStr = ''
    for key, value in segDepth.items():
        resStr += key+'\t'+str(value)+'\n'
    segFile = open('{}_seg.txt'.format(args.sampleName), 'w')
    segFile.write(resStr)
    segFile.close()

if __name__ == "__main__":
    main()