from MutationMatrixGenerator import MutationMatrixGenerator
import numpy as np
import pandas as pd
import pickle
from pathlib import Path

class MutationTableGenerator(MutationMatrixGenerator):
    def __init__(self, species,  possible_species, exon_file, gene_annotations,  output_dir, trinuc_prob_file,N,Nb):
        super().__init__(species,  possible_species, exon_file, gene_annotations, output_dir)
        self.output_syn_dir = Path(f'{output_dir}/syn_target')
        self.output_nonsyn_dir = Path(f'{output_dir}/nonsyn_target')
        self.trinuc_prob_df = pd.read_csv(trinuc_prob_file, index_col=0, sep='\t')
        self.r_WS = (2*Nb)/(1-np.exp(-2*Nb))  #fixation probability
        self.r_SW = N*(1 - np.exp(2*Nb/N))/(1 - np.exp(2*Nb))

    def find_gc_bias(self,ref, alt):
    #need to define mapping between the bases and multiply to add the gc bias 
        if ref in ['G','C']:
            ref_bias = 1 #1 represents strong
        else:
            ref_bias = 0 #0 represents weak
        if alt in ['G','C']:
            alt_bias = 1
        else:
            alt_bias = 0
        if ref_bias == 1 and alt_bias == 0:
            gc_bias = self.r_SW
        elif ref_bias == 0 and alt_bias == 1:
            gc_bias = self.r_WS
        else:
            gc_bias = 1
        return gc_bias 

  

    #generate and output all non syn mutations and save to one file
    #seperatly generate all syn muts 
    def generate_mutations(self, seq, trinucs, type):
        results=[]
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon not in self.codon_table:
                continue
            for m in range(3):
                seq_pos = i + m
                if seq_pos >= len(trinucs):
                    continue
                ref_base = seq[seq_pos]
                for alt in self.bases:
                    if alt == ref_base:
                        continue
                    alt_codon = list(codon)
                    alt_codon[m] = alt
                    alt_codon = ''.join(alt_codon)
                    mut_type = self.find_mutation_type(codon, alt_codon)
                    if mut_type == type:
                        context_index = trinucs[seq_pos]
                        if context_index != -1:
                            trinuc_string = self.trinucleotides[context_index]
                            trinuc_prob = self.trinuc_prob_df.at[trinuc_string, alt]
                            gc_prob = self.find_gc_bias(ref_base,alt)
                            prob = gc_prob*trinuc_prob
                            results.append([seq_pos, trinuc_string, alt, prob])
        
        results_df = pd.DataFrame({
            'gene_pos': [x[0] for x in results],
            'trinuc_string': [x[1] for x in results],
            'alt': [x[2] for x in results],
            'prob': [x[3] for x in results]   
        })
        
        return results_df



    #overwritting to save table of muts + porbs, rather than matrix
    def run(self):
        x=0
        self.load_data()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.output_syn_dir.mkdir(parents=True, exist_ok=True)
        self.output_nonsyn_dir.mkdir(parents=True, exist_ok=True)
        print(f'number of genes in filtered df {len(self.all_exons_df['gene'].unique())}')
        with open(self.output_dir / 'processed_gene_seqs.pkl', "ab") as f:
            for gene, gene_df in self.all_exons_df.groupby("gene"):
                #save start coordinate also? 
                seq, trinucs = self.extract_gene_seq_trinucs(gene_df)
                record=(gene,seq,trinucs)
                pickle.dump(record, f, protocol=pickle.HIGHEST_PROTOCOL)
                if len(seq) % 3 == 0:
                    results_nonsyn = self.generate_mutations(seq, trinucs, 1) #1 is for non syn
                    if not results_nonsyn.empty:
                        results_nonsyn.to_csv(self.output_nonsyn_dir / f"{gene}",index=False)
                    results_syn = self.generate_mutations(seq, trinucs, 0) #0 is for  syn
                    if not results_syn.empty:
                        results_syn.to_csv(self.output_syn_dir / f"{gene}",index=False)
                else:
                    x+=1
                    print(x)