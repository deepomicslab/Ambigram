import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--del', dest='de', required=True, help='')
parser.add_argument('-i', '--inv', dest='inv', required=True, help='')
parser.add_argument('-p', '--dup', dest='dup', required=True, help='')
parser.add_argument('-t', '--trans', dest='trans', required=True, help='')
parser.add_argument('-m', '--map', dest='map', required=True, help='')
parser.add_argument('-s', '--seg', dest='seg', required=True, help='')
parser.add_argument('-H', '--hap', dest='hap', required=True, help='')
args = parser.parse_args()

del_df = pd.read_table(args.de)
inv_df = pd.read_table(args.inv)
dup_df = pd.read_table(args.dup)
trans_df = pd.read_table(args.trans)
bps_map = pd.read_table(args.map, names=['chrom', 'before', 'after'])
segs = pd.read_table(args.seg)
with open(args.hap) as fin:
    h1, h2 = [line[:-1].split() for line in fin.readlines()]

for i in del_df.index:
	after = bps_map.loc[lambda df: df.before >= del_df.loc[i, 'Start'] + 28460000 - 1 - 5].loc[lambda df: df.before <= del_df.loc[i, 'Start'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		del_df.at[i, 'Start'] = after.iloc[0]
	else:
		print('Fail')
	after = bps_map.loc[lambda df: df.before >= del_df.loc[i, 'End'] + 28460000 - 1 - 5].loc[lambda df: df.before <= del_df.loc[i, 'End'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		del_df.at[i, 'End'] = after.iloc[0]
	else:
		print('Fail')

for i in inv_df.index:
	after = bps_map.loc[lambda df: df.before >= inv_df.loc[i, 'Start'] + 28460000 - 1 - 5].loc[lambda df: df.before <= inv_df.loc[i, 'Start'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		inv_df.at[i, 'Start'] = after.iloc[0]
	else:
		print('Fail')
	after = bps_map.loc[lambda df: df.before >= inv_df.loc[i, 'End'] + 28460000 - 1 - 5].loc[lambda df: df.before <= inv_df.loc[i, 'End'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		inv_df.at[i, 'End'] = after.iloc[0]
	else:
		print('Fail')


for i in dup_df.index:
	after = bps_map.loc[lambda df: df.before >= dup_df.loc[i, 'Start'] + 28460000 - 1 - 5].loc[lambda df: df.before <= dup_df.loc[i, 'Start'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		dup_df.at[i, 'Start'] = after.iloc[0]
	else:
		print('Fail')
	after = bps_map.loc[lambda df: df.before >= dup_df.loc[i, 'End'] + 28460000 - 1 - 5].loc[lambda df: df.before <= dup_df.loc[i, 'End'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		dup_df.at[i, 'End'] = after.iloc[0]
	else:
		print('Fail')

for i in trans_df.index:
	after = bps_map.loc[lambda df: df.before >= trans_df.loc[i, 'StartA'] + 28460000 - 1 - 5].loc[lambda df: df.before <= trans_df.loc[i, 'StartA'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		trans_df.at[i, 'StartA'] = after.iloc[0]
	else:
		print('Fail')
	after = bps_map.loc[lambda df: df.before >= trans_df.loc[i, 'EndA'] + 28460000 - 1 - 5].loc[lambda df: df.before <= trans_df.loc[i, 'EndA'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		trans_df.at[i, 'EndA'] = after.iloc[0]
	else:
		print('Fail')
	after = bps_map.loc[lambda df: df.before >= trans_df.loc[i, 'StartB'] + 28460000 - 1 - 5].loc[lambda df: df.before <= trans_df.loc[i, 'StartB'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		trans_df.at[i, 'StartB'] = after.iloc[0]
	else:
		print('Fail')
	after = bps_map.loc[lambda df: df.before >= trans_df.loc[i, 'EndB'] + 28460000 - 1 - 5].loc[lambda df: df.before <= trans_df.loc[i, 'EndB'] + 28460000 - 1 + 5].after
	if len(after.unique()) == 1:
		trans_df.at[i, 'EndB'] = after.iloc[0]
	else:
		print('Fail')
    

th1 = []
th2 = []
for row in segs.itertuples():
	del_entry = del_df.loc[lambda df: df.Start == row.start].loc[lambda df: df.End == row.end]
	inv_entry = inv_df.loc[lambda df: df.Start == row.start].loc[lambda df: df.End == row.end]
	dup_entry = dup_df.loc[lambda df: df.Start == row.start].loc[lambda df: df.End == row.end]
	trans_entryA = trans_df.loc[lambda df: df.StartA == row.start].loc[lambda df: df.EndA == row.end]
	trans_entryB = trans_df.loc[lambda df: df.StartB == row.start].loc[lambda df: df.EndB == row.end]
	if len(del_entry) > 0:
		if del_entry.Chr.iloc[0] == 'mhc_1':
			th2.append(str(row.ID) + '+')
		elif del_entry.Chr.iloc[0] == 'mhc_2':
			th1.append(str(row.ID) + '+')
	elif len(inv_entry) > 0:
		if inv_entry.Chr.iloc[0] == 'mhc_1':
			th1.append(str(row.ID) + '-')
			th2.append(str(row.ID) + '+')
		elif inv_entry.Chr.iloc[0] == 'mhc_2':
			th1.append(str(row.ID) + '+')
			th2.append(str(row.ID) + '-')
	elif len(dup_entry) > 0:
		if dup_entry.Chr.iloc[0] == 'mhc_1':
			th1.extend([str(row.ID) + '+'] * dup_entry.Duplications.iloc[0])
			th2.append(str(row.ID) + '+')
		elif dup_entry.Chr.iloc[0] == 'mhc_2':
			th1.append(str(row.ID) + '+')
			th2.extend([str(row.ID) + '+'] * dup_entry.Duplications.iloc[0])
	elif len(trans_entryA) > 0:
		if trans_entryA.ChrB.iloc[0] == 'mhc_1':
			th1.append(str(segs.loc[lambda df: df.start == trans_entryA.StartB.iloc[0]].ID.iloc[0]) + '+')
			th2.append(str(row.ID) + '+')
		elif trans_entryA.ChrB.iloc[0] == 'mhc_2':
			th1.append(str(row.ID) + '+')
			th2.append(str(segs.loc[lambda df: df.start == trans_entryA.StartB.iloc[0]].ID.iloc[0]) + '+')
	elif len(trans_entryB) > 0:
		if trans_entryB.ChrA.iloc[0] == 'mhc_1':
			th1.append(str(segs.loc[lambda df: df.start == trans_entryB.StartA.iloc[0]].ID.iloc[0]) + '+')
			th2.append(str(row.ID) + '+')
		elif trans_entryB.ChrA.iloc[0] == 'mhc_2':
			th1.append(str(row.ID) + '+')
			th2.append(str(segs.loc[lambda df: df.start == trans_entryB.StartA.iloc[0]].ID.iloc[0]) + '+')
	else:
		th1.append(str(row.ID) + '+')
		th2.append(str(row.ID) + '+')

print(' '.join(th1))
print(' '.join(th2))

n_del_correct = 0
for row in del_df.itertuples():
    try:
        seg_id = segs.loc[lambda df: df.start == row.Start].iloc[0].ID
        v = str(seg_id) + '+'
        if (h1 + h2).count(v) == 1:
            n_del_correct += 1
    except:
        pass

n_inv_correct = 0
for row in inv_df.itertuples():
    try:
        seg_id = segs.loc[lambda df: df.start == row.Start].iloc[0].ID
        v = str(seg_id) + '-'
        if v in h1:
            idx = h1.index(v)
            a = [int(i[:-1]) for i in h1[(idx - 1):(idx + 2)]]
        else:
            idx = h2.index(v)
            a = [int(i[:-1]) for i in h2[(idx - 1):(idx + 2)]]
        if a[1] - a[0] == 1 and a[2] - a[1] == 1:
            n_inv_correct += 1
    except:
        pass

n_dup_correct = 0
for row in dup_df.itertuples():
    try:
        seg_id = segs.loc[lambda df: df.start == row.Start].iloc[0].ID
        v = str(seg_id) + '+'
        if (h1 + h2).count(v) - 1 == row.Duplications:
            n_dup_correct += 1
    except:
        pass

n_trans_correct = 0
for row in trans_df.itertuples():
    try:
        segA_id = segs.loc[lambda df: df.start == row.StartA].iloc[0].ID
        segB_id = segs.loc[lambda df: df.start == row.StartB].iloc[0].ID
        vA = str(segA_id) + '+'
        vB = str(segB_id) + '+'
        h = h1 + h2
        th = th1 + th2
        vA_idx = [idx for idx, e in enumerate(h) if e == vA]
        vB_idx = [idx for idx, e in enumerate(h) if e == vB]
        vA_tidx = [idx for idx, e in enumerate(th) if e == vA]
        vB_tidx = [idx for idx, e in enumerate(th) if e == vB]
        for i in vA_tidx:
            for j in vA_idx:
                if th[(i - 1):(i + 2)] == h[(j - 1):(j + 2)]:
                    n_trans_correct += 0.5
        for i in vB_tidx:
            for j in vB_idx:
                if th[(i - 1):(i + 2)] == h[(j - 1):(j + 2)]:
                    n_trans_correct += 0.5
    except:
        pass
print(f'Del: {n_del_correct}')
print(f'Inv: {n_inv_correct}')
print(f'Dup: {n_dup_correct}')
print(f'Trans: {n_trans_correct}')