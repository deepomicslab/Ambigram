import sys
import pandas as pd

sample_name = sys.argv[1]

with open(f'{sample_name}/haps/results') as fin:
	th1, th2 = [line[:-1].split() for line in fin.readlines()[:2]]

with open(f'{sample_name}/haps/{sample_name}.haps') as fin:
	h1, h2 = [line[:-1].split() for line in fin.readlines()]

df = []

th = th1 + th2
th_s = [int(i[:-1]) for i in th]
th_su = sorted(set(th_s))
h = h1 + h2
h_s = [int(i[:-1]) for i in h]
h_su = sorted(set(h_s))

for s in th_su:
	t_count = th_s.count(s)
	count = h_s.count(s)
	if t_count > 2:
		df.append((s, 'dup', t_count, count, count - t_count))
	elif t_count == 2:
		df.append((s, 'norm', t_count, count, count - t_count))
	else:
		df.append((s, 'del', t_count, count, count - t_count))

df = pd.DataFrame(df, columns=['seg', 'type', 't_count', 'count', 'diff'])

df.to_csv(f'{sample_name}/haps/compare', index=False, sep='\t')

norm_c = {'more': 0, "less": 0}
dup_c = {'more': 0, "less": 0}
del_c = {'more': 0, "less": 0}
for row in df.itertuples():
    if row.type == 'norm':
        if row.diff > 0:
            norm_c['more'] += 1
        elif row.diff < 0:
            norm_c['less'] += 1
    if row.type == 'dup':
        if row.diff > 0:
            dup_c['more'] += 1
        elif row.diff < 0:
            dup_c['less'] += 1
    if row.type == 'del':
        if row.diff > 0:
            del_c['more'] += 1
        elif row.diff < 0:
            del_c['less'] += 1
diff = pd.DataFrame({'norm' : norm_c, 'dup': dup_c, 'del': del_c}).T
diff.to_csv(f'{sample_name}/haps/compare.diff', index=True, sep='\t')