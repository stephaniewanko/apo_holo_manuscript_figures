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
from figure_functions import *	

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

order_all[order_all.s2ang_x < 0] = 0
order_all[order_all.s2ang_y < 0] = 0
order_all[order_all.s2calc_x < 0] = 0


#MERGE
merged_order_all = merge_apo_holo_df(order_all)
print(merged_order_all.head())

#All Order Parameters by Residue Type
merged_order_all['s2calc_x'] = merged_order_all['s2calc_x'].clip(lower=0)
merged_order_all['s2calc_y'] = merged_order_all['s2calc_y'].clip(lower=0) 

#SUBSET
merged_order_all_polar = merged_order_all[merged_order_all['resn_x'].isin(['R','N','D','C','Q','E', 'H', 'K', 'S', 'T', 'Y'])]
merged_order_all_nonpolar = merged_order_all[merged_order_all['resn_x'].isin(['V','W','P','F','M','L','I','G','A'])]

#All Order Parameter Distribution Plots
make_dist_plot_AH(merged_order_all['s2calc_x'], merged_order_all['s2calc_y'], 's2calc', 'Number of Residues', 'Apo v. Holo s2calc (Entire Protein)', 'AH_s2calc_entireprotein')
make_dist_plot_AH(merged_order_all['s2ortho_x'], merged_order_all['s2ortho_y'], 's2ortho', 'Number of Residues', 'Apo v. Holo s2ortho (Entire Protein)', 'AH_s2ortho_entireprotein')

#STATS
print('Difference of s2calc on only Polar Side Chains between Holo/Apo [Entire Protein]')
paired_ttest(merged_order_all_polar['s2calc_x'], merged_order_all_polar['s2calc_x'])

print('Difference of s2calc on only nonpolar Side Chains between Holo/Apo [Entire Protein]')
paired_ttest(merged_order_all_nonpolar['s2calc_x'], merged_order_all_nonpolar['s2calc_x'])


#FIGURE 
fig = plt.figure()
f, axes = plt.subplots(1, 4, sharey=True, sharex=True)

p1 = sns.boxenplot(merged_order_all_polar['s2calc_x'], orient='v', 
ax=axes[0]).set(xlabel='Polar Bound', ylabel='S2calc Order Parameter')
p2 = sns.boxenplot(merged_order_all_polar['s2calc_y'], orient='v', ax=axes[1]).set(xlabel='Polar Unbound', ylabel='')
p3 = sns.boxenplot(merged_order_all_nonpolar['s2calc_x'], orient='v', ax=axes[2]).set(xlabel='NonPolar Bound', ylabel='')
p4 = sns.boxenplot(merged_order_all_nonpolar['s2calc_y'], orient='v', ax=axes[3]).set(xlabel='NonPolar Unbound', ylabel='')
fig.savefig('FullProtein_OP_byResidueType.png')
