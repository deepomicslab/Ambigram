def parse_ilp_result(filename):
    d = {}
    with open(filename) as fin:
        for line in fin.readlines()[1:]:
            line = line[:-1].split()
#             print(line)
            dd = {int(line[1][1:]) : float(line[2])}
#             print(dd)
            d.update(dd)
    return d

# def parse_ilp_result(filename):
#     d = {}
#     with open(filename) as fin:
#         for line in fin.readlines():
#             if re.search('x[0-9]+', line):
#                 line = list(filter(None, re.split('\s+| |\*+|=+', line)))
#                 print(f"{line[1][1:]} : {line[2]}")
#                 d.update({int(line[1][1:]) : float(line[2])})
#     return d

def generate_balanced_lh(outfilename, infilename, d):
    with open(outfilename, 'w') as fout:
        with open(infilename) as fin:
#             n_seg = 0
            junc_id = 0
            for line in fin.readlines():
                line = line[:-1].split()
                if line[0] == 'SEG':
                    seg_id = int(line[1].split(':')[1]) - 1
                    try:
                        seg_cp = d[seg_id]
                    except KeyError:
#                         print(seg_id)
                        seg_cp = 0
#                     print(f"{seg_id} : {seg_cp}")
                    junc_id += 1
                    fout.write(' '.join(line[:-2]) + f' {seg_cp} {line[-1]}\n')
#                     print(' '.join(line[:-1]) + f' {seg_cp} {line[-1]}\n')
                elif line[0] == 'JUNC':
                    try:
                        junc_cp = d[junc_id]
                    except KeyError:
#                         print(junc_id)
                        junc_cp = 0
#                     print(f"{junc_id} : {junc_cp}")
                    junc_id += 1
                    fout.write(' '.join(line[:-3]) + f' {junc_cp} ' + ' '.join(line[-2:]) + '\n')
#                     print(' '.join(line[:-2]) + f' {junc_cp} ' + ' '.join(line[-2:]) + '\n')
                else:
                    fout.write(' '.join(line) + '\n')
#                     print(' '.join(line) + '\n')