import pandas as pd
from IntronMatrix import IntronMatrix
import os
'''
chunk and read introns (removing cpg and dups as go
Then save to matrixs with targets and muts)'''

#introns from cactus allignments


class IntronMatrixNoCpG(IntronMatrix):
    def __init__(self,seq_file, output_dir,species, target_species, possible_species):
        super().__init__(seq_file, output_dir,species, target_species, possible_species)
        
    def run(self):
        os.makedirs(self.output_dir,exist_ok=True)
        header = ['chr', 'seq_start','seq_end', 'target_seq','seq', 'chr_copy','start', 'end', 'gene','gene_name','strand' ]
        chunk_iter = pd.read_csv(self.seq_file, chunksize=1000, sep='\t', names=header)
        
        for chunk in chunk_iter:
            chunk.drop(columns=['chr_copy','gene','gene_name'],inplace=True)
            chunk.dropna(inplace=True)
            chunk['start'] = chunk['start'].astype(int)
            chunk['seq_start'] = chunk['seq_start'].astype(int)
            chunk['end'] = chunk['end'].astype(int)

            #getting error of overlap too small and seq = Na

            chunk['seq'] = chunk.apply(lambda row: row['seq'][row['start']-row['seq_start']: row['end']-row['seq_start']],axis=1)
            chunk['target_seq']=chunk.apply(lambda row: row['target_seq'][row['start']-row['seq_start']: row['end']-row['seq_start']],axis=1)
            
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
        return
        

