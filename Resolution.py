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

make_dist_plot_AH(AH_pairs['Holo_Res'], AH_pairs['Apo_Res'], 'Resolution', 'Number of Structures', 'Apo v. Holo s2calc (Entire Protein)', 'FullResolution')
print('Difference of s2calc on only Polar Side Chains between Holo/Apo [Entire Protein]')
paired_ttest(AH_pairs['Apo_Res'], AH_pairs['Holo_Res'])
