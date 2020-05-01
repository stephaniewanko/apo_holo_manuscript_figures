
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

#READ IN FILES
os.chdir('/Users/stephaniewankowicz/Dropbox/Fraser_Rotation/qfit_output/200225/')
path=os.getcwd()
all_files = glob.glob(path + "/*qFit_qfit_RMSF.csv")
li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=0, header=0)
    li.append(df)
RMSF = pd.concat(li, axis=0, ignore_index=True)
RMSF['PDB'] = RMSF['PDB_name'].astype(str).str[0:4]
print(len(all_files))

RMSF_merged = merge_apo_holo_df(RMSF)

make_dist_plot_AH(RMSF_merged['rmsf_x'], RMSF_merged['rmsf_y'],'RMSF Entire Protein', 'Number of Residues', 'RMSF across entrie protein', 'RMSF')
print('Difference of RMSF on Side Chains between Holo/Apo [Entire Protein]')
paired_ttest(RMSF_merged['rmsf_x'], RMSF_merged['rmsf_y'])


RMSF_summary = pd.DataFrame()
n = 1
for i in RMSF['PDB'].unique():
    tmp = RMSF[RMSF['PDB'] == i]
    RMSF_summary.loc[n, 'PDB'] = i
    RMSF_summary.loc[n, 'Num_Residues'] = len(tmp.index)
    RMSF_summary.loc[n, 'Num_Alt_Loc'] = len(tmp[tmp['rmsf']>0].index)
    if tmp.rmsf.ge(0).any() == True:
        RMSF_summary.loc[n, 'Alt_Loc'] = 1
    else:
        RMSF_summary.loc[n, 'Alt_Loc'] = 0
    RMSF_summary.loc[n, 'Apo_Holo'] = tmp['Apo_Holo'].unique()
    RMSF_summary.loc[n, 'Average_RMSF'] = tmp['rmsf'].mean()
    n += 1

RMSF_summary['per_altloc'] = RMSF_summary['Num_Alt_Loc'] / RMSF_summary['Num_Residues']
RMSF_summary['Num_Single'] = RMSF_summary['Num_Residues'] - RMSF_summary['Num_Alt_Loc']


RMSF_summary_m= merge_apo_holo_df(RMSF_summary)

RMSF_summary_m['Apo_Holo_Multi_Diff'] = RMSF_summary_m['per_altloc_x'] - RMSF_summary_m['per_altloc_y']

same = len((RMSF_summary_m[RMSF_summary_m['Apo_Holo_Multi_Diff']==0]).index)
gain = len((RMSF_summary_m[RMSF_summary_m['Apo_Holo_Multi_Diff']>0]).index)
loss = len((RMSF_summary_m[RMSF_summary_m['Apo_Holo_Multi_Diff']<0]).index)

#ENTIRE PROTEIN
x = range(3)
print(same)
print(gain)
print(loss)
p1 = plt.bar(x[0], same, width=0.7)
p2 = plt.bar(x[1], gain, width=0.7)
p3 = plt.bar(x[2], loss, width=0.7)

#plt.title('Holo - Apo Pairs')
plt.ylabel('Number of Pairs (Holo-Apo)')
plt.xticks(x, ('Same','Gain', 'Loss'))
plt.show()

x = range(3)
f, axes = plt.subplots(1, 3, sharey=True, sharex=True)
p1 = sns.boxenplot(RMSF_summary['per_altloc'], orient='v', ax=axes[0]).set(xlabel='All', ylabel='% Residues with Alt Loc')
p2 = sns.boxenplot(RMSF_summary[RMSF_summary['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[1]).set(xlabel='Unbound', ylabel='')
p3 = sns.boxenplot(RMSF_summary[RMSF_summary['Apo_Holo']=='Holo']['per_altloc'], orient='v', ax=axes[2]).set(xlabel='Bound', ylabel='')
plt.show()


os.chdir('/Users/stephaniewankowicz/Dropbox/Fraser_Rotation/qfit_output/200225/')
path=os.getcwd()

all_files = glob.glob(path + "/*5.0_rmsf_subset.csv")

li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=0, header=0)
    li.append(df)

close_RMSF = pd.concat(li, axis=0, ignore_index=True)
print(len(all_files))


merged_close_RMSF = merge_apo_holo_df(close_RMSF)

merged_close_sum_RMSF['Percent_Holo_Close'] = merged_close_sum_RMSF['Num_Alt_Loc_x']/ merged_close_sum_RMSF['Num_Residues_x']
merged_close_sum_RMSF['Percent_Apo_Close'] = merged_close_sum_RMSF['Num_Alt_Loc_y']/ merged_close_sum_RMSF['Num_Residues_y']


merged_close_sum_RMSF['Apo_Holo_Multi_Diff'] = merged_close_sum_RMSF['Percent_Holo_Close'] - merged_close_sum_RMSF['Percent_Apo_Close']


same = len((merged_close_sum_RMSF[merged_close_sum_RMSF['Apo_Holo_Multi_Diff']==0]).index)
gain = len(merged_close_sum_RMSF[merged_close_sum_RMSF['Apo_Holo_Multi_Diff']>0].index)
loss = len((merged_close_sum_RMSF[merged_close_sum_RMSF['Apo_Holo_Multi_Diff']<0]).index)

x = range(3)
print(same)
print(gain)
print(loss)
p1 = plt.bar(x[0], same, width=0.7)
p2 = plt.bar(x[1], gain, width=0.7)
p3 = plt.bar(x[2], loss, width=0.7)

plt.ylabel('Number of Pairs')
plt.xticks(x, ('Same','Gain', 'Loss'))
plt.show()