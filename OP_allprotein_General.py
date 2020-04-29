#packages
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from scipy import stats
import matplotlib.patches as mpatches

#reference files
os.chdir('/Users/stephaniewankowicz/Dropbox/Fraser_Rotation/qfit_output/')
pairs = pd.read_csv('ligand_supplementary_table1.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

os.chdir('/Users/stephaniewankowicz/Dropbox/Fraser_Rotation/qfit_output/')
AH_key=pd.read_csv('qfit_AH_key_191218.csv')

#read in files
os.chdir('/Users/stephaniewankowicz/Dropbox/Fraser_Rotation/qfit_output/200225/')
path=os.getcwd()


all_files = glob.glob(path + "/*qFit_methyl.out")

li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    df['PDB'] = filename[69:73]
    li.append(df)

order_all = pd.concat(li, axis=0, ignore_index=True)
print(len(all_files))


#All Order Parameters by Residue Type
order_all = order_all.merge(AH_key, on=['PDB'])
order_all_holo = order_all[order_all['Apo_Holo']=='Holo']
order_all_apo = order_all[order_all['Apo_Holo']=='Apo']

test = order_all.merge(AH_pairs, left_on='PDB', right_on='Holo')
merged_order_all = test.merge(order_all_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi']) 
merged_order_all = merged_order_all.drop_duplicates()

merged_order_all['s2calc_x'] = merged_order_all['s2calc_x'].clip(lower=0)
merged_order_all['s2calc_y'] = merged_order_all['s2calc_y'].clip(lower=0) 

#SUBSET
merged_order_all_polar = merged_order_all[merged_order_all['resn_x'].isin(['R','N','D','C','Q','E', 'H', 'K', 'S', 'T', 'Y'])]

merged_order_all_nonpolar = merged_order_all[merged_order_all['resn_x'].isin(['V','W','P','F','M','L','I','G','A'])]


#STATS
print('Difference of s2calc on only Polar Side Chains between Holo/Apo [Entire Protein]')
print(stats.ttest_rel(merged_order_all_polar['s2calc_x'], merged_order_all_polar['s2calc_y']))

print('Holo Mean:')
print(merged_order_all_polar['s2calc_x'].mean())

print('Apo Mean:')
print(merged_order_all_polar['s2calc_y'].mean())

print('Holo Median:')
print(merged_order_all_polar['s2calc_x'].median())

print('Apo Mean:')
print(merged_order_all_polar['s2calc_y'].median())


print('Difference of s2calc on only nonpolar Side Chains between Holo/Apo [Entire Protein]')
print(stats.ttest_rel(merged_order_all_nonpolar['s2calc_x'], merged_order_all_nonpolar['s2calc_y']))

print('Holo Mean:')
print(merged_order_all_nonpolar['s2calc_x'].mean())

print('Apo Mean:')
print(merged_order_all_nonpolar['s2calc_y'].mean())

print('Holo Median:')
print(merged_order_all_nonpolar['s2calc_x'].median())

print('Apo Mean:')
print(merged_order_all_nonpolar['s2calc_y'].median())


#FIGURE
fig = plt.figure()
f, axes = plt.subplots(1, 4, sharey=True, sharex=True)

p1 = sns.boxenplot(merged_order_all_polar['s2calc_x'], orient='v', 
ax=axes[0]).set(xlabel='Polar Bound', ylabel='Scalc Order Parameter')
p2 = sns.boxenplot(merged_order_all_polar['s2calc_y'], orient='v', ax=axes[1]).set(xlabel='Polar Unbound', ylabel='')
p3 = sns.boxenplot(merged_order_all_nonpolar['s2calc_x'], orient='v', ax=axes[2]).set(xlabel='NonPolar Bound', ylabel='')
p4 = sns.boxenplot(merged_order_all_nonpolar['s2calc_y'], orient='v', ax=axes[3]).set(xlabel='NonPolar Unbound', ylabel='')
plt.show()
figure.savefig('FullProtein_OP_byResidueType.png')
