import numpy as np
import pandas as pd
import pickle
import glob
import os
import copy

'''
inputs
syn and non syn directories with files labelled by gene with pos, prob (including trinuc and gc bias) and alt
start speci file- pickle df containing gene, filtered cds sequence in direction of transcription and tirnuc contexts
cbase output-estimates of lambda n and lambda s 
PIPELINE
multinomial resitrbute total mutations across genes using cbase outputs
within each gene redistirbute syn and non syn muts seperatly using inputted probs, and multinomial sampling
apply muts to seqs in dicts, and chnage correponding tirnucs 
OUTPUT
mutated df with gene, cds seq and trinucs
'''


class NeutralEvolutionSimulator:
    def __init__(self, n_subs, start_speci_file, cbase_output, syn_probs_dir, non_syn_probs_dir, output,gene_name_map):
        self.Nsubs = n_subs
        self.start_speci_file = start_speci_file
        self.cbase_output = cbase_output
        self.syn_probs_dir = syn_probs_dir
        self.non_syn_probs_dir = non_syn_probs_dir
        self.output = output
        #self.gene_name_map = pd.read_csv(gene_name_map,index_col=1)

        # Constants
        self.bases = ['A', 'C', 'G', 'T']
        self.trinucleotides = [a + b + c for a in self.bases for b in self.bases for c in self.bases]
       

        

    def load_pickle_records_to_df(self,path):
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

        return pd.DataFrame(rows,columns=['gene','seq','trinucs']) 

    def load_cbase_data(self):
        cbase_output_df = pd.read_csv(self.cbase_output, index_col=0)
        #map from gene name to ensembl ids
        #cbase_output_df = cbase_output_df.join(self.gene_name_map, how='inner').set_index('gene')
        
        #filter to only genes in reocnstruction 
        # Get gene file *basenames* (strip directories)
        genes_syn = [os.path.basename(p) for p in glob.glob(f'{self.syn_probs_dir}/*')]
        genes_nonsyn = [os.path.basename(p) for p in glob.glob(f'{self.non_syn_probs_dir}/*')]

        genes = set(genes_syn).union(set(genes_nonsyn))
        print(f'number of valid recon genes is {len(genes)}')
        self.cbase_output_df = cbase_output_df[cbase_output_df.index.isin(genes)]
        print(cbase_output_df[~cbase_output_df.index.isin(genes)])

    #need to check if this keeps correct order
    def distribute_muts_across_genes(self):
        probs = self.cbase_output_df[['lambda_s', 'lambda_n']].to_numpy().flatten()
        total = probs.sum()
        genes_probs= (probs / total).tolist()
        list_of_genes = self.cbase_output_df.index.values
        print(f'number of genes with cbase output is {len(list_of_genes)}')
        mutations = np.random.multinomial(self.Nsubs, genes_probs)
        #from index 0 in steps of 2
        muts_synon = mutations[::2]
        #from index 1 in steps of 2
        muts_non_synon = mutations[1::2]
        data = {'gene': list_of_genes, 'synon_muts': muts_synon, 'non_synon_muts': muts_non_synon}
        df = pd.DataFrame(data).set_index('gene')
        
        #keeping only mutated genes
        df = df[(df['synon_muts']+df['non_synon_muts'])>0]
        return df

    #be careful that the reults are in the same order as the exon file, 
    # as we then also append direct entries from exon file
    def apply_trinucs_map(self,trinuc_index,alt):
        trinuc= self.trinucleotides[trinuc_index]
        alt_trinuc = trinuc[0] + alt + trinuc[2]
        return self.trinucleotides.index(alt_trinuc)


    def apply_muts_per_gene(self, sim_df, gene, gene_df):
        gene_seq = gene_df['seq'].iloc[0]
        gene_trinucs = gene_df['trinucs'].iloc[0]

        alt_map = dict(zip(sim_df['gene_pos'], sim_df['alt']))
        seq_list = list(gene_seq)
        mutated_seq = ''.join(alt_map.get(i, base) for i, base in enumerate(gene_seq))
        '''
        for idx, trinuc in enumerate(gene_trinucs):
            if idx in sim_df['gene_pos'].values:
                gene_trinucs[idx] = self.apply_trinucs_map(trinuc,alt_map.get(idx) )
        '''
        mut_trinucs = copy.deepcopy(gene_trinucs)
        for pos,alt in alt_map.items():
            trinuc = gene_trinucs[pos]
            mut_trinucs[pos] = self.apply_trinucs_map(trinuc,alt)
        
        record=(gene,mutated_seq,mut_trinucs)
        
        with open(self.output, "ab") as f:
            pickle.dump(record, f, protocol=pickle.HIGHEST_PROTOCOL)
        return

    def run(self):
        #clear file if already exists, as well append
        open(self.output, 'w').close()

        self.load_cbase_data()
        muts_per_gene_df = self.distribute_muts_across_genes()
        mutated_genes = muts_per_gene_df.index.unique()
        print(f'number of mutated genes in simulation is {len(mutated_genes)}')
        #load gene seqs, to simulate muts
        start_speci_df =  self.load_pickle_records_to_df(self.start_speci_file)
        mutated_start_speci_df =  start_speci_df[start_speci_df['gene'].isin(mutated_genes)]
        
        for gene, gene_row in mutated_start_speci_df.groupby('gene'):
            syn_gene_df = pd.read_csv(f'{self.syn_probs_dir}/{gene}')
            non_syn_gene_df = pd.read_csv(f'{self.non_syn_probs_dir}/{gene}')
            #check as previously eccidently saved empty dfs
            if syn_gene_df.empty or non_syn_gene_df.empty:
                continue
            syn_gene_df = syn_gene_df.copy()
            non_syn_gene_df = non_syn_gene_df.copy()

            Ngene_syn = muts_per_gene_df.loc[gene,'synon_muts']
            Ngene_non_syn = muts_per_gene_df.loc[gene,'non_synon_muts'] 
            

            syn_sum = syn_gene_df['prob'].sum()
            syn_gene_df['norm_prob'] = syn_gene_df['prob']/syn_sum
            syn_gene_df['sim_mut'] = np.random.multinomial(Ngene_syn, syn_gene_df['norm_prob'].values)
            
            non_syn_sum = non_syn_gene_df['prob'].sum()
            non_syn_gene_df['norm_prob'] = non_syn_gene_df['prob']/non_syn_sum
            non_syn_gene_df['sim_mut'] = np.random.multinomial(Ngene_non_syn, non_syn_gene_df['norm_prob'].values)

            gene_simulated_subs_df = pd.concat([syn_gene_df, non_syn_gene_df]).sort_index()
            gene_no0_simulated_subs_df = gene_simulated_subs_df[gene_simulated_subs_df.sim_mut != 0]

            if not gene_no0_simulated_subs_df.empty:
                self.apply_muts_per_gene(gene_no0_simulated_subs_df, gene, gene_row)

        # Save non-mutated genes
        unmutated_start_speci_df =  start_speci_df[~start_speci_df['gene'].isin(mutated_genes)]
        with open(self.output, "ab") as f:
            for idx, row in unmutated_start_speci_df.iterrows():   # list of dicts per row
                rec = (row['gene'], row['seq'], row['trinucs'])
                pickle.dump(rec, f, protocol=pickle.HIGHEST_PROTOCOL)
            



