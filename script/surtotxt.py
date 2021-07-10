import sys

in_sur = sys.argv[1]
out_txt = sys.argv[2]

in_lh = open(in_sur)
f = open(out_txt)
for line in f.readlines():
    a = line.split("\t")
    #     in_lh.write(JUNC H:16:+ H:11:+ 969.5 -1 U B)
    in_lh.write("JUNC H:" + a[0] + ":+ H:" + a[0] + ":+ 969.5 -1 U B\n")
    in_lh.write("JUNC H:" + a[0] + ":+ H:" + a[0] + ":- 969.5 -1 U B\n")
    in_lh.write("JUNC H:" + a[0] + ":- H:" + a[0] + ":- 969.5 -1 U B\n")