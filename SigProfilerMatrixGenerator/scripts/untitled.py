







trans = mut_count_all['6144'][mut_count_all['6144'].index.str.contains('T:')]
untrans = mut_count_all['6144'][mut_count_all['6144'].index.str.contains('U:')]
bitrans = mut_count_all['6144'][mut_count_all['6144'].index.str.contains('B:')]
nontrans = mut_count_all['6144'][mut_count_all['6144'].index.str.contains('N:')]

a = bitrans//2
b = bitrans%2
trans = trans.append(a)
trans  = trans.groupby(trans.index.str[2:]).sum()
trans = trans.rename(index=lambda s: "T:"+s)
untrans = untrans.append(a)
untrans = untrans.groupby(untrans.index.str[2:]).sum()
untrans = untrans.rename(index=lambda s:"U:"+s)

mut_count_all['4608'] = pd.concat([finalCollapse, nontrans])
for x in list(b.index):
	mutAbrrev = x.split(":")[1]
	tempMat = pd.concat([trans.loc["T:"+mutAbrrev],untrans.loc["U:"+mutAbrrev]], axis=1)
	for z,y in zip(tempMat.idxmax(axis=1), list(tempMat.index)):
		tempMat.loc[y,z] += b.loc[x,y]

	finalCollapse = finalCollapse.append(tempMat, ignore_index=True)








biTransNonTrans =  mut_count_all['6144'][mut_count_all['6144'].index.str.contains('B:') | mut_count_all['6144'].index.str.contains('N:')]
biTransNonTrans = biTransNonTrans.groupby(biTransNonTrans.index.str[2:]).sum()
biTransNonTrans = biTransNonTrans.rename(index=lambda s: "N:"+s)
mut_count_all['4608'] = pd.concat([transUntrans, biTransNonTrans])









