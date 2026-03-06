import pandas as pd
from collections import Counter
import numpy as np
import os

'''
INPUT
intron file,
OUTPUT
intron trinuc mut matrix, sub matrix and trinuc count matrix
PIPELINE
iterate through genome 
reverse complement 
save muts and correpsonding trinucs
save counts of each trinuc in introns
 '''




class IntronMatrix:
    def __init__(self, seq_file, output_dir, species, target_species, possible_species):
        self.bases = ['A', 'C', 'G', 'T']
        self.trinucleotides = [a + b + c for a in self.bases for b in self.bases for c in self.bases]
        self.rc_table = str.maketrans({'A':'T', 'C':'G', 'G':'C', 'T':'A'})
        self.seq_file = seq_file
        self.possible_species = possible_species
        self.species=species
        self.target_species = target_species
        
        self.output_dir = output_dir
        self.mut_df =  pd.DataFrame(0.0, index=self.trinucleotides, columns=self.bases)
        self.trinuc_df = pd.DataFrame(0, index=self.trinucleotides, columns=['count'])
        #load dfs for muts and trnucs here 

    def reverse_complement(self, seq):
        bases = list(seq)
        comp_bases=[]
        for base in bases:
            if base in self.complement.keys():
                comp_base = self.complement[base]
            else:
                comp_base = base
            comp_bases.append(comp_base)
        comp_seq = ''.join(comp_bases)
        return comp_seq[::-1]

    def save_trinucs_bin(self, row_seq):
        counts = Counter(tri for i in range(len(row_seq) - 2) if (tri := row_seq[i:i+3]) in self.trinucleotides)
        
        # Convert to DataFrame with trinucleotide strings as the index
        df = pd.DataFrame.from_dict(counts, orient='index', columns=['count'])
        # Make sure all trinucleotides are represented, filling missing ones with 0
        df = df.reindex(self.trinucleotides, fill_value=0)
        #print('check count right number of trinucs')
        #print(sum(1 for c in seq if c in self.bases))
        #print(df.sum())
        # Efficient assignment: assign the entire column for one bin - check order CHECX
        
        self.trinuc_df['count'] += df['count']
        
        return 

    def find_mutations(self, human_seq, HCLCA_seq):
        l = min(len(human_seq), len(HCLCA_seq))
        for k in range(1,l-1):
            if human_seq[k] != HCLCA_seq[k]:
                trinuc = HCLCA_seq[k-1:k+2]
                alt = human_seq[k]
                if trinuc in self.trinucleotides and alt in self.bases:
                    self.mut_df.loc[trinuc,alt] +=1
        return

    def run(self):
        #os.mkdir(self.output_dir)
        header = ['chr', 'start', 'end'] + self.possible_species + [ 'gene', 'strand']
        chunk_iter = pd.read_csv(self.seq_file, chunksize=1000, sep='\t', names=header)
        
        for chunk in chunk_iter:
            chunk = chunk[[self.species, self.target_species, 'strand']]
            chunk.rename(columns={self.species:'seq', self.target_species:'target_seq'},inplace=True)
            
            chunk['seq'] = chunk.apply(lambda row: row['seq'].translate(self.rc_table)[::-1]
                 if row['strand'] == '-' else row['seq'],axis=1)
            chunk['target_seq'] = chunk.apply(lambda row: 
                row['target_seq'].translate(self.rc_table)[::-1] if row['strand'] == '-' else row['target_seq'],axis=1)
            chunk.apply(lambda row: self.save_trinucs_bin(row['seq']),axis=1)
            chunk.apply(lambda row: self.find_mutations(row['target_seq'], row['seq']),axis=1)
        
        if (self.mut_df.index==self.trinuc_df.index).all():
            normalised_df = self.mut_df.div(self.trinuc_df['count'],axis=0)
            normalised_df.to_csv(f'{self.output_dir}/rates', sep='\t')
        self.mut_df.to_csv(f'{self.output_dir}/muts', sep='\t')
        self.trinuc_df.to_csv(f'{self.output_dir}/trinucs', sep='\t')
        
