#Figure Functions
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def paired_ttest(holo_col, apo_col):
	print(stats.ttest_rel(holo_col, apo_col))

	print('Holo Mean:')
	print(holo_col.mean())

	print('Apo Mean:')
	print(apo_col.mean())

	print('Holo Median:')
	print(holo_col.median())

	print('Apo Mean:')
	print(apo_col.median())


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
