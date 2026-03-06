import pickle
import pandas as pd
from collections import Counter
target='hg38'
anc='Anc4'
episode = f'{target}_{anc}'
#input parameters
seq_file = f'/home/maria/run_simulations_cactus_clean/{episode}/sim_auxiliary/{anc}_processed_gene_seqs.pkl'
output_file=f'/home/maria/run_simulations_cactus_clean/{episode}/sim_auxiliary/{anc}_trinucs.pkl'

bases = ['A', 'C', 'G', 'T']
trinucleotides = [a + b + c for a in bases for b in bases for c in bases]
trinuc_df = pd.DataFrame(0, index=trinucleotides, columns=['count'])

def load_pickle_records_to_df(path):
    rows = []
    with open(path, "rb") as f:
        while True:
            try:
                obj = pickle.load(f)
            except EOFError:
                break

            if isinstance(obj, list):     # came from chunked writer
                rows.extend(obj)
            else:                          # came from per-record writer
                rows.append(obj)

    return pd.DataFrame(rows,columns=['gene','seq','trinucs']).set_index('gene') 

df=load_pickle_records_to_df(seq_file)
gene_trinuc_dict={}
#already save syn and non syn targets per gene in simulation runs, should check how will it matches cbase occs.pkl file for the synonymous
for gene,row in df.iterrows():
    #eg [12,13,21,22], count instance of each 
    gene_trinuc_df=trinuc_df.copy()
    trinucs_counts=Counter(row['trinucs'])
    for i in range(64):
        if i in trinucs_counts:
            gene_trinuc_df.loc[trinucleotides[i],'count']=trinucs_counts[i]
    
    gene_trinuc_dict[gene]=gene_trinuc_df
    #update

print(gene_trinuc_dict)
with open(output_file, "ab") as f:
    pickle.dump(gene_trinuc_dict, f, protocol=pickle.HIGHEST_PROTOCOL)