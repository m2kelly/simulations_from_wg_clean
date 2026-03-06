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
cbase_syn_target_dfs=pickle.load(open(f'/home/maria/run_simulations_cactus_clean/scripts/2_run_cbase/hg38_Anc4/hg38_Anc4/Anc4->hg38/genes/occs_per_gene.pkl', 'rb'))

syn_files=glob.glob(f'{syn_target_dir}/*')
for gene_file in syn_files[:1]:
    target_df=pd.read_csv(gene_file,index_col=0)
    gene=gene_file.split('/')[-1]
    cbase_df=cbase_syn_target_dfs[gene].transpose()
    print(cbase_df)
    '''
    target df has colns trinuc string alt and prob, want to output df with indexes trinuc_string
    columns alt, and entries the count of occurences in target_df for each trinuc_stirng,alt combo
    
    '''
    
    out = (
        target_df
        .groupby(['trinuc_string','alt'])
        .size()
        .unstack(fill_value=0)
        .reindex(index=CONTEXTS, columns=BASES, fill_value=0)
    )
    print(out)