import pickle
import pandas as pd
from collections import Counter
import glob
target='hg38'
anc='Anc4'
episode = f'{target}_{anc}'
BASES = ['A','C','G','T']
CONTEXTS = [a+b+c for a in BASES for b in BASES for c in BASES]

#cbase outputs occurences, which should be possible syn muts, compare to my outputted targets
syn_target_dir='/home/maria/run_simulations_cactus_clean/hg38_Anc4/sim_auxiliary/syn_target'
nonsyn_target_dir='/home/maria/run_simulations_cactus_clean/hg38_Anc4/sim_auxiliary/nonsyn_target'
qval_file='/home/maria/run_simulations_cactus_clean/scripts/2_run_cbase/hg38_Anc4/hg38_Anc4/Anc4->hg38/q_values_hg38_Anc4.txt'
output_file='/home/maria/run_simulations_cactus_clean/hg38_Anc4/sim_auxiliary/Anc4_syn_frac.pkl'
#only run for cbase signficant genes
qval_df=pd.read_csv(qval_file,skiprows=1,sep='\t',index_col=0)
signf_genes=qval_df[qval_df['q_phi_neg']<0.05].index.tolist()

gene_dicts={}
for gene in signf_genes:
    try:
        syn_df=pd.read_csv(f'{syn_target_dir}/{gene}',index_col=0)
    except:
        continue
    try: 
        nonsyn_df=pd.read_csv(f'{nonsyn_target_dir}/{gene}',index_col=0)
    except:
        continue
    
    '''
    target df has colns trinuc string alt and prob, want to output df with indexes trinuc_string
    columns alt, and entries the count of occurences in target_df for each trinuc_stirng,alt combo
    
    '''
    
    syn = (
        syn_df
        .groupby(['trinuc_string','alt'])
        .size()
        .unstack(fill_value=0)
        .reindex(index=CONTEXTS, columns=BASES, fill_value=0)
    )
    nonsyn = (
        nonsyn_df
        .groupby(['trinuc_string','alt'])
        .size()
        .unstack(fill_value=0)
        .reindex(index=CONTEXTS, columns=BASES, fill_value=0)
    )
    total=syn+nonsyn
    syn_frac=syn.div(total)
    syn_frac.fillna(value=0.0,inplace=True)
    gene_dicts[gene]=syn_frac


with open(output_file,"ab") as f:
    pickle.dump(gene_dicts, f, protocol=pickle.HIGHEST_PROTOCOL)

  