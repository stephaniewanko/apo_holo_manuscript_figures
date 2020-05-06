#Figure Functions
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd

#reference files
os.chdir('/Users/stephaniewankowicz/Dropbox/Fraser_Rotation/qfit_output/')
pairs = pd.read_csv('ligand_supplementary_table1.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

os.chdir('/Users/stephaniewankowicz/Dropbox/Fraser_Rotation/qfit_output/')
AH_key=pd.read_csv('qfit_AH_key_191218.csv')

def paired_ttest(holo_col, apo_col):
	print(stats.ttest_rel(holo_col, apo_col))

	print('Holo Mean:')
	print(holo_col.mean())

	print('Apo Mean:')
	print(apo_col.mean())

	print('Holo Median:')
	print(holo_col.median())

	print('Apo Median:')
	print(apo_col.median())



def ind_ttest(col_1, col_2):
	print(stats.ttest_ind(col_1, col_2))

	print('Column 1 Mean:')
	print(col_1.mean())

	print('Column 2 Mean:')
	print(col_2.mean())

	print('Column 1 Median:')
	print(col_1.median())

	print('Column 2 Median:')
	print(col_2.median())



def merge_apo_holo_df(df):
	df = df.merge(AH_key, on = ['PDB'])
	df_holo = df[df['Apo_Holo'] == 'Holo']
	df_apo = df[df['Apo_Holo'] == 'Apo']

	test = df_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
	df_merged = test.merge(df_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
	df_merged = df_merged.drop_duplicates()
	return df_merged



def make_dist_plot_AH(holo_col, apo_col, x_label, y_label, title, out_name):
    fig = plt.figure()
    sns.distplot(holo_col, kde = False, label = 'Holo')
    sns.distplot(apo_col, kde = False, label = 'Apo')
    plt.xlabel(x_label)
    plt.legend()
    plt.ylabel(y_label)
    plt.title(title)
    fig.savefig(out_name + '.png')



def make_boxenplot_chem(low_col, high_col, xlabel_low, xlabel_high, ylabel, out_name):
	fig = plt.figure()
	x = range(2)
	f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
	p1 = sns.boxenplot(low_col, orient='v', ax=axes[0]).set(xlabel = xlabel_low, ylabel = ylabel)
	p2 = sns.boxenplot(high_col, orient='v', ax=axes[1]).set(xlabel = xlabel_high, ylabel = '')
	plt.savefig(out_name + '.png')
	#plt.show()

