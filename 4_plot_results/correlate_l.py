import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import os
from sklearn.linear_model import LinearRegression
plt.rcParams.update({'font.size': 14})


ANC='Anc4'
TARGET='hg38'


'''
correlate changes in l with s het and dn ds
'''

#global variables
episode = f'{TARGET}_{ANC}'
cancer_genes= '/home/maria/run_simulations/data/Cosmic_CancerGeneCensus_names.txt'
gene_ids = '/home/maria/cactus_target_size/auxillary/gene_name_id.csv'
WEIGHTS_FILE = '/home/maria/data/signature_weights/calc_mut_weights.csv'
DIR = f'/home/maria/run_simulations_cactus_clean/{episode}'
TARGET_FILE= f'{DIR}/output/{TARGET}'
ANC_FILE= f'{DIR}/output/{ANC}'
PLOT_DIR= f'{DIR}/plots'
shet_file='/home/maria/papers/shet/Supplementary_Table_1.txt'
cancer_type='BRCA'
qvals_file=f'/home/maria/data/cbase_MC3/q_values_MC3_original_{cancer_type}.txt'

qvals_file=f'/home/maria/run_simulations_cactus_clean/scripts/2_run_cbase/{episode}/{episode}/{ANC}->{TARGET}/q_values_{episode}.txt'

#make directory
if not os.path.exists(PLOT_DIR):
    os.makedirs(PLOT_DIR)

#OPEN DFS TO BE SUED MULTIPLE TIMES
GENE_ID_DF = pd.read_csv(gene_ids)
FILTERED_SUM = pd.read_csv(WEIGHTS_FILE, index_col=0, names=['norm_weight'])
cancer_genes_df = pd.read_csv(cancer_genes, names=['gene_name'])
#convert cancer genes to ensembl gene ids
cancer_id_df = cancer_genes_df.merge(GENE_ID_DF, on=['gene_name'], how='inner')
CANCER_GENES_LIST = cancer_id_df['gene'].tolist()
shet_df = pd.read_csv(shet_file,sep='\t',index_col=0)  #use 's_het_drift'
qval_df = pd.read_csv(qvals_file, sep='\t',skiprows=1,index_col=0)
#to compare to shet
#GENE_ID_DF.set_index(keys='gene_name',inplace=True)
#qval_df = qval_df.join(GENE_ID_DF,how='inner')
#qval_df.set_index(keys='gene',inplace=True)


#correlate qvals with delta l in human, regress on number of observed muts
sig = 'SBS5'
anc_df = pd.read_csv(f'{ANC_FILE}/l_{sig}',index_col=0,names=['l_anc'])
target_df = pd.read_csv(f'{TARGET_FILE}/l_{sig}',index_col=0,names=['l_target'])
l_df = anc_df.join(target_df,how='inner')
l_df['l_diff'] = (l_df['l_anc']-l_df['l_target'])/(l_df['l_anc']+l_df['l_target'])

qval_df['subs'] = qval_df['k_obs']+qval_df['s_obs']+ qval_df['m_obs']

q_l_df = l_df.join(qval_df, how='inner')
#restrict to signficant q value genes
q_l_df= q_l_df[q_l_df['q_phi_neg']<0.05]
print(len(q_l_df))

# Multilinear regression: l_diff ~ d(m+k)/ds with offset subs
X = q_l_df[['d(m+k)/ds', 'subs']].values
#X = q_l_df['d(m+k)/ds'].values.reshape(-1,1)
y = q_l_df['l_diff'].values

model = LinearRegression()
model.fit(X, y)

print(f"Multilinear Regression Results:")
print(f"Coefficient d(m+k)/ds: {model.coef_[0]}")
print(f"Coefficient subs: {model.coef_[1]}")
print(f"Intercept: {model.intercept_}")
print(f"R-squared: {model.score(X, y)}")

#if diff>0 -> reduced susceptibility
#small d(m+k)/ds -> strong germline selection

#add a regression of mean simulated suscepitbility - target susceptibility 
#vs dn/ds 
#plot simulations
no_of_sims=100
sim_file = f"{DIR}/output/sim_{{x}}/l_{sig}"
dfs=[]
for i in range(no_of_sims):
    l_df = pd.read_csv(sim_file.format(x=i), header=None, index_col=0)
    l_df.columns = [i]
    dfs.append(l_df)
sim_df = pd.concat(dfs,axis=1)
sim_df['l_sim'] = sim_df.mean(axis=1)

sim_q_df = sim_df.join(q_l_df,how='inner')
print(stats.spearmanr(sim_q_df['l_sim']-sim_q_df['l_target'],sim_q_df['d(m+k)/ds']))
