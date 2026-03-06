import pandas as pd
from collections import Counter
import numpy as np
import pickle
from pathlib import Path
from joblib import Parallel, delayed
from concurrent.futures import ProcessPoolExecutor
import sys


'''
INPUT df 
from 1_ints_bins_sub_cpg
with refernce genome and cpg island, both orginal coordinates intersected 
+df with cpg islands subtracted
+intersect chr bin positions? 
OUTPUT picked dictionarys per chromsome, key is bin and entry is df of trinuc mutations
+ df of trinuc occurances
PROCESS
extract the cpg islands per bin, 
then bin genome by chromosome
do not include chrX or chrY
for each sub extrct trinuc island + count trinucs per bin 

OUTPUT
dictionary, keys are chromsomes
entries are dfs
indexes are poss muts eg TTT->A , length 192 , do not allow TTT->T 
columns are bins, starting with 0
'''


class MatrixSub:
    def __init__(self, output, bin_seq_file,name):
        self.name=name
        self.bin_seq_file = bin_seq_file
        self.output=f'{output}/{name}'
        self.bases = ['A', 'C', 'G', 'T']
        self.c_bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        self.trinucleotides = [a + b + c for a in self.bases for b in self.bases for c in self.bases]
        self.all_possible_muts = [f"{trinuc}->{base}" for trinuc in self.trinucleotides for base in self.bases if base!=trinuc[1]]
        Path(self.output).mkdir(parents=True, exist_ok=True)
    
    #input se
    def process_row(self, bin_start, seq, target_seq):
        for i in range(1,len(seq)-1): #cannot do edge bases, as do not have trinucs 
            ref = seq[i]
            alt = target_seq[i]
            if ref != alt: #sub has happened, save
                trinuc = seq[i-1:i+2]
                if (trinuc in self.trinucleotides) and (alt in self.bases):
                    self.mut_df.loc[f"{trinuc}->{alt}", bin_start] +=1
    
        return 


    def save_trinucs_bin(self,bin_start, row_seq):
        allowed = self.trinucleotides
        s = row_seq
        counts = Counter()
        n = len(s) - 2
        for i in range(n):
            tri = s[i:i+3]
            if tri in allowed:
                counts[tri] += 1

        counts = Counter(tri for i in range(len(row_seq) - 2) if (tri := row_seq[i:i+3]) in self.trinucleotides)
        
        #format to be have index AAA->C, AAA->G, AAA->T, rather than just AAA for each trinuc

        trinuc_counts = {mut: counts.get(mut[:3], 0) for mut in self.all_possible_muts}
        self.trinuc_df.loc[self.all_possible_muts, bin_start] += pd.Series(trinuc_counts)
        
        
        return 

    def process_bins(self,bin_name,group, allowed):
        counts = Counter()
        for s in group["seq"].tolist():
            n = len(s) - 2
            for i in range(n):
                tri = s[i:i+3]
                if tri in allowed:
                    counts[tri] += 1
        # one write per bin
        self.trinuc_df.loc[list(counts.keys()), bin_name] += pd.Series(counts)
        
        #process subs
        group.apply(lambda row: self.process_row(row['bin_name'], row['seq'],row['target_seq']),axis=1)
        return



    def mutations_structure(self,chr_df):
        no_bins = max(chr_df['bin_name']) +1
        starts = list(range(no_bins))
        self.mut_df = pd.DataFrame(0, index=self.all_possible_muts, columns=starts)
        self.trinuc_df = pd.DataFrame(0, index=self.trinucleotides, columns=starts)
        allowed = self.trinucleotides  # set

        #when paralise crashes memory and does not update parent self.trinuc_df
        #proess by bin group, subs and count trinucs
        #Parallel(n_jobs=multiprocessing.cpu_count())(delayed(self.process_bins)((name, group, allowed) for name, group in chr_df.groupby("bin_name", sort=False)))
        #with ProcessPoolExecutor(max_workers=4) as executor:
            #executor.map(self.process_bins, ((name, group, allowed) for name, group in chr_df.groupby("bin_name", sort=False)))

        for name, group in chr_df.groupby("bin_name", sort=False):
            self.process_bins(name, group, allowed)

        #convert trinucs df to muts df indices
        tri_for_mut = [m[:3] for m in self.all_possible_muts]  # e.g. 'AAA' for 'AAA->C'
        muts_vals = self.trinuc_df.reindex(tri_for_mut).to_numpy(copy=False)
        self.trinuc_muts_df = pd.DataFrame(
            muts_vals,
            index=self.all_possible_muts,
            columns=self.trinuc_df.columns
        )
        return
    
    def load_bins(self):
        self.bins_df = pd.read_csv(self.bin_seq_file, sep='\t', names = ['chr', 'seq_start','seq_end', 'target_seq','seq', 'chr_copy1','bin_start', 'bin_end', 'chr_copy2', 'start', 'end'])
        #test
        #self.bins_df = self.bins_df.iloc[:1000]
        

        self.bins_df.drop(columns=['chr_copy1', 'chr_copy2','bin_end','seq_end'],inplace=True)
        self.bins_df['start'] = self.bins_df['start'].astype(int)
        self.bins_df['seq_start'] = self.bins_df['seq_start'].astype(int)
        self.bins_df['end'] = self.bins_df['end'].astype(int)
        self.bins_df = self.bins_df[self.bins_df['chr'].isin([f'chr{i}' for i in range(1,23)])]
        return
    
    #notusing
    def load_chr(self,chr):
        self.bins_df = pd.read_csv(self.bin_seq_file, sep='\t', names = ['chr', 'seq_start','seq_end', 'target_seq','seq', 'chr_copy1','bin_start', 'bin_end', 'chr_copy2', 'start', 'end'])
        self.bin_df = self.bins_df[self.bins_df['chr']==chr]
        
        self.bins_df.drop(columns=['chr_copy1', 'chr_copy2','bin_end','seq_end'],inplace=True)
        self.bins_df['start'] = self.bins_df['start'].astype(int)
        self.bins_df['seq_start'] = self.bins_df['seq_start'].astype(int)
        self.bins_df['end'] = self.bins_df['end'].astype(int)
        
        return
    
    

    def run(self):
        self.load_bins()
        self.muts_dict = {}
        self.trinuc_dict = {}
        #self.bins_df.sort_values(by = ['start'], inplace=True)
        for chr, chr_df in self.bins_df.groupby('chr'):
            
            print(f'processing {chr}')
            #chr_df['seq'] = chr_df.apply(lambda row: row['seq'][row['start']-row['seq_start']: row['end']-row['seq_start']],axis=1)
            #chr_df['target_seq']=chr_df.apply(lambda row: row['target_seq'][row['start']-row['seq_start']: row['end']-row['seq_start']],axis=1)
            #remove seqs of length 1
            chr_df = chr_df[chr_df['end'] - chr_df['start'] >3]
            #cut seqns
            starts = (chr_df['start'].to_numpy() - chr_df['seq_start'].to_numpy()).astype(int)
            ends   = (chr_df['end'].to_numpy()   - chr_df['seq_start'].to_numpy()).astype(int)

            seqs = chr_df['seq'].to_list()
            tseqs = chr_df['target_seq'].to_list()

            chr_df = chr_df.copy()
            try:
                chr_df['seq'] = [s[a:b] for s, a, b in zip(seqs, starts, ends)]
                chr_df['target_seq'] = [s[a:b] for s, a, b in zip(tseqs, starts, ends)]
            except Exception as e:
                print(f"Error processing {chr}: {e}")
                for i, s in enumerate(seqs):
                    if not isinstance(s, str):
                        print(f"Index {i}: {s} ({type(s)})")
                
                sys.exit("Condition not met — exiting.")
            
            chr_df['bin_name'] = chr_df['bin_start'] // 100000
            #chr_df['bin_name'] = chr_df['bin_start'].apply(lambda x :0 if x==0 else int(str(x)[:-5])) #remove 10^5
            self.mutations_structure(chr_df)
            
            self.muts_dict[chr] = self.mut_df
            self.trinuc_dict[chr] = self.trinuc_muts_df
            
            
        with open(f'{self.output}/whole_genome_dict_{self.name}_unsmoothed_muts.pkl', 'wb') as fp:
            pickle.dump(self.muts_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f'{self.output}/whole_genome_dict_{self.name}_unsmoothed_occs.pkl', 'wb') as fp:
            pickle.dump(self.trinuc_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
        


