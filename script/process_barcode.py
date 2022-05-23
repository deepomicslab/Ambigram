import argparse

def readSeg(segDir):
    segFile = open(segDir, "r")
    segs = []
    for line in segFile.readlines():
        info = line.split("\t")[0]
        chr = info.split(":")[0]
        pos = info.split(":")[1]
        segs.append([chr, int(pos.split("-")[0]), int(pos.split("-")[1]), len(segs)+1])
    return segs

def readBarcode(barcodeDir, segs):
    barcodeFile = open(barcodeDir, "r")
    group =  [[] for i in range(len(segs))]
    for line in barcodeFile.readlines():
        info = line.strip('\n').split("\t")
        chr = info[0]
        if chr[0] != 'c':
            chr = "chr"+chr
        pos1 = int(info[1])
        pos2 = int(info[2])
        code = info[3]
        start = -1
        end = -1
        min1 = float('inf')
        min2 = float('inf')
        for i in range(0,len(segs)):            
            seg = segs[i]
            if chr!=seg[0]:
                continue
            if i==0 and pos1 <= seg[1]:
                start = i
            elif i==len(segs)-1 and pos2 >= seg[2]:
                end = i
            else:
                if abs(seg[1]-pos1)<min1:
                    start = i
                    min1 = abs(seg[1]-pos1)
                if abs(seg[2]-pos2)<min2:
                    end = i
                    min2 = abs(seg[2]-pos2)
        if start>end or start not in range(0,len(segs)) or end not in range(0,len(segs)):
            continue
        for i in range(start,end+1):
            group[i].append(code)
    for arr in group:
        arr = list(set(arr))
    return group

def getLinkWeight(group,i,j):
    if i>=j:
        return 0
    res = set(group[i])
    for idx in range(i+1,j+1):
        res = res & set(group[idx])
    return len(res)

def barcode2juncs(segDir, barcodeDir, juncDir):
    # segDir = "./segments/COLO829_1_seg.txt"
    segs = readSeg(segDir)
    print(segs)
    # barcodeDir = "./third_expr/10x/coverage/10x_30x.bed"
    group = readBarcode(barcodeDir,segs)
    # print(group)

    intervals = []
    source = 0
    for i in range(1,len(segs)):
        if segs[i][0]!=segs[source][0]:
            intervals.append([source,i-1])
            source = i
    if(source<len(segs)):
        intervals.append([source,len(segs)-1])
    print(intervals)
    links = []
    # link segments
    for interval in intervals:
        for i in range(interval[0],interval[1]):
            for j in range(i+1,interval[1]+1):
                w = getLinkWeight(group,i,j)*(j-i)
                # print('seg{}-seg{}:{}'.format(i+1,j+1,w))
                links.append([i+1,j+1,w])
    links.sort(key=lambda x: x[2], reverse=True)
    for info in links:
        print('seg{}-seg{}:{}'.format(info[0],info[1],info[2]))
    res = ""
    for i in range(0,5):
        for seg in range(links[i][0],links[i][1]):
            res += str(seg)+"+ "
        res += str(links[i][1])+"+\n"
    juncFile = open(juncDir, "w")
    juncFile.write(res)
    # translocations
    # trans = []
    # for i in range(0,len(intervals)):
    #     for j in range(i+1,len(intervals)):
    #         for u in range(intervals[i][0],intervals[i][1]+1):
    #             for v in range(intervals[j][0],intervals[j][1]+1):
    #                 w = getLinkWeight(group,u,v)
    #                 trans.append([u+1,v+1,w])
    # trans.sort(key=lambda x: x[2])
    # print("Possible translocations:")
    # for info in trans:
    #     print('seg{}-seg{}:{}'.format(info[0],info[1],info[2]))   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate .lh file by integrating sv.txt and seg.txt.')
    parser.add_argument('-bed', '--bed_file', dest='bedPath', required=True, help='Path to BED file from BarcodeExtract')
    parser.add_argument('-seg', '--seg_file', dest='segPath', required=True, help='Path to SEG file')
    parser.add_argument('-s', '--sample_name', dest='sampleName', required=False, default = 'sample', help='Sample name for output file')
    args = parser.parse_args()
    juncDir = '{}.juncs'.format(args.sample)
    barcode2juncs(args.bedPath, args.segPath, juncDir)