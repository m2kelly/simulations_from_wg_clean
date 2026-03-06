import numpy as np
from scipy.stats import stats
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

cancer_genes = '/home/maria/run_simulations/data/Cosmic_CancerGeneCensus_names.txt'
#cancer_genes='/home/maria/data/Census_TSG.txt'
gene_ids = '/home/maria/cactus_target_size/auxillary/gene_name_id.csv'
#homo_file = '/home/maria/cactus_target_size/output/hg38'
#HCLCA_file = '/home/maria/cactus_target_size/output/anc'
#removing cpgs and dups
anc='Anc1'
target='Anc3'
homo_file = f'/home/maria/run_simulations_cactus_clean/{target}_{anc}/output_no_div_M_n/{target}'
HCLCA_file=f'/home/maria/run_simulations_cactus_clean/{target}_{anc}/output_no_div_M_n/{anc}'

homo_file = f'/home/maria/run_simulations_cactus_clean/{target}_{anc}/output_x_100000/{target}'
HCLCA_file=f'/home/maria/run_simulations_cactus_clean/{target}_{anc}/output_x_100000/{anc}'

homo_file = f'/home/maria/run_simulations_cactus_clean/{target}_{anc}/output/{target}'
HCLCA_file=f'/home/maria/run_simulations_cactus_clean/{target}_{anc}/output/{anc}'

mu_file = '/home/maria/synon_mut_rates/auxillary_nocpg_nodup/neutral_exon_muts_norm_5/muts_per_genes.tsv'

sig = 'SBS5'
l_n_df1 = pd.read_csv(f'{homo_file}/l_{sig}',index_col=0,names=['l_n1'])
l_n_df2 = pd.read_csv(f'{HCLCA_file}/l_{sig}',index_col=0,names=['l_n2'])
l_n_df = l_n_df1.join(l_n_df2,how='inner')
#l_n_df['l_n_diff'] = (l_n_df['l_n2'] - l_n_df['l_n1'])/(l_n_df['l_n2'] + l_n_df['l_n1'])
l_n_df['l_n_diff'] = (l_n_df['l_n2'] - l_n_df['l_n1'])

mu_df = pd.read_csv(mu_file,index_col=0)
combined_df = l_n_df.join(mu_df,how='inner')
combined_df.dropna(inplace=True)

print(f'l_n*mu reduction = {100*(combined_df['l_n_diff']*combined_df['mean_mu']).sum()/(combined_df['l_n2']*combined_df['mean_mu']).sum()} %')
print(f'l_n reduction = {100*l_n_df['l_n_diff'].sum()/(l_n_df['l_n2']).sum()} %')


#restricting to cancer genes
gene_id_df = pd.read_csv(gene_ids)
cancer_genes_df = pd.read_csv(cancer_genes, names=['gene_name'])
#convert cancer genes to ensembl gene ids
cancer_id_df = cancer_genes_df.merge(gene_id_df, on=['gene_name'], how='inner')
cancer_genes_list = cancer_id_df['gene'].tolist()

#cancer_combined_df = combined_df
cancer_combined_df = combined_df[combined_df.index.isin(cancer_genes_list)]
cancer_l_n_df = l_n_df[l_n_df.index.isin(cancer_genes_list)]
print(f'cancer l_n reduction = {100*cancer_l_n_df['l_n_diff'].sum()/(cancer_l_n_df['l_n2']).sum()} %')
print(f'cancer l_n*mu reduction = {100*(cancer_combined_df['l_n_diff']*cancer_combined_df['mean_mu']).sum()/(cancer_combined_df['l_n2']*cancer_combined_df['mean_mu']).sum()} %')

