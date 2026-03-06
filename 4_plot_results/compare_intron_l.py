import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import os
from sklearn.linear_model import LinearRegression
plt.rcParams.update({'font.size': 14})


ANC='fullTreeAnc213'
TARGET='fullTreeAnc212'


'''
correlate changes in l with s het and dn ds
'''

#global variables
episode = f'{TARGET}_{ANC}'
#cancer_genes= '/home/maria/run_simulations/data/Cosmic_CancerGeneCensus_names.txt'
#gene_ids = '/home/maria/cactus_target_size/auxillary/gene_name_id.csv'
#WEIGHTS_FILE = '/home/maria/data/signature_weights/calc_mut_weights.csv'
DIR = f'/home/maria/run_simulations_cactus_clean/{episode}'
TARGET_FILE= f'{DIR}/output_intron_matrix/{TARGET}/l_intron'
ANC_FILE= f'{DIR}/output_intron_matrix/{ANC}/l_intron'
output_plot=f'{DIR}/plots/l_intron.png'

#correlate qvals with delta l in human, regress on number of observed muts

anc_df = pd.read_csv(ANC_FILE,index_col=0,names=['l_anc'])
target_df = pd.read_csv(TARGET_FILE,index_col=0,names=['l_target'])
l_df = anc_df.join(target_df,how='inner')
l_df['l_diff'] = (l_df['l_anc']-l_df['l_target'])/(l_df['l_anc']+l_df['l_target'])
print(len(l_df[l_df['l_diff']>0]))
print(len(l_df[l_df['l_diff']<0]))
print(len(l_df))
l_df=l_df[l_df['l_diff']!=0]
plt.hist(l_df['l_diff'],bins=100)
plt.xlabel(r'$\Delta l_n$, intron matrix per branch')
plt.ylabel('frequency')
plt.title(episode)
plt.axvline(x=0,label='no change')
mean=l_df['l_diff'].mean()
plt.axvline(x=mean,label=f'mean change {mean:3f}',color='red')
plt.legend()
plt.savefig(output_plot)
plt.show()