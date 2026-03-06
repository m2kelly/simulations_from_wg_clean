import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import Counter
target='Anc3'
anc = 'Anc1'
episode=f'{target}_{anc}'
sim_bedfile=f'/home/maria/run_simulations_cactus_clean/{episode}/output/simulated_genes_{{x}}.pkl'
'''
load seq pickle and count no of gc bases/gc ending codons
+no of syn and non syn muts
'''
def count_total_gc_content_in_gene(seq):
    counts=Counter(seq)
        
    return counts['A'], counts['C'], counts['G'], counts['T']
    
def count_oppurtunities(seq):
    M_n = 0
    M_s =0 
    AT_ending =0
    GC_ending = 0
    syn_AT_ending =0
    syn_GC_ending = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon not in codon_table:
            continue
        if codon[0:2] in ['AC', 'GT','GC', 'GG', 'TC', 'CT', 'CC', 'CG']: #first 2 bases of all 4 fold degenerate sites
            if codon[2] in ['A','T']:
                syn_AT_ending += 1
            if codon[2] in ['G','C']:
                syn_GC_ending +=1
        #checking translational bias
        if codon[2] in ['A','T']:
            AT_ending += 1
        if codon[2] in ['G','C']:
            GC_ending +=1

        #calculating number of synonymous and non synonymous oppurutnities
        for m in range(3):
            seq_pos = i + m
            #just to be sure , should be fine with i range
            if seq_pos >= len(seq):
                continue
            ref_base = seq[seq_pos]
            
            for alt in self.bases:
                
                if alt == ref_base:
                    continue
                alt_codon = list(codon)
                alt_codon[m] = alt
                alt_codon = ''.join(alt_codon)
                mut_type = self.find_mutation_type(codon, alt_codon)
                if mut_type == 1:
                    M_n +=1
                elif mut_type == 0:
                    M_s += 1

    return M_n, M_s, AT_ending, GC_ending, syn_AT_ending, syn_GC_ending


output_folder = '/home/maria/cactus_target_size/output_primates_extended'
analysis_file = '/home/maria/cactus_target_size/primate_extended/auxillaryprimate_genome_heuristics'
possible_species = ['Homo_sapiens','fullTreeAnc105','fullTreeAnc106','fullTreeAnc107','fullTreeAnc108','fullTreeAnc109']



plt.rcParams.update({'font.size': 18})

bases = ['A', 'C', 'G', 'T']

#extracting list of cancer gene ids 
cancer_genes = '/home/maria/run_simulations/data/Cosmic_CancerGeneCensus_names.txt'
gene_ids = '/home/maria/cactus_target_size/auxillary/gene_name_id.csv'
gene_id_df = pd.read_csv(gene_ids)
cancer_genes_df = pd.read_csv(cancer_genes, names=['gene_name'])
#convert cancer genes to ensembl gene ids
cancer_id_df = cancer_genes_df.merge(gene_id_df, on=['gene_name'], how='inner')
cancer_genes_list = cancer_id_df['gene'].tolist()


speci_gc_percentage ={}
speci_non_syn_percentage = {}
speci_GC_ending_percentage = {}
for speci in possible_species:
    speci_df = pd.read_csv(f'{analysis_file}/{speci}', index_col =0)
    #restricting to cancer genes lists
    speci_df = speci_df[speci_df.index.isin(cancer_genes_list)]
    speci_sum_df = speci_df.sum()
    speci_gc_percentage[speci] = (speci_sum_df['C'] + speci_sum_df['G'])/(speci_sum_df['A'] + speci_sum_df['T']+speci_sum_df['C'] + speci_sum_df['G'])
    speci_non_syn_percentage[speci] = (speci_sum_df['M_n'])/(speci_sum_df['M_n'] + speci_sum_df['M_s'])
    speci_GC_ending_percentage[speci] = (speci_sum_df['syn_GC_ending'])/(speci_sum_df['syn_AT_ending'] + speci_sum_df['syn_GC_ending'])
    '''
    y = np.array([speci_sum_df[base] for base in bases])
    plt.figure()
    plt.pie(y, labels = bases)
    plt.title(speci)
    '''
    #plt.show() 


print(speci_gc_percentage)
point_sizes = {speci:30+i*30 for i,speci in enumerate(possible_species)}

species_dict ={'Homo_sapiens':'Homo_sapiens','fullTreeAnc105':'Anc1','fullTreeAnc106':'Anc2','fullTreeAnc107':'Anc3','fullTreeAnc108':'Anc4','fullTreeAnc109':'Anc5'}

plt.figure()
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(species_dict[speci], speci_gc_percentage[speci], label=species_dict[speci], s = size)
plt.title('Cancer Genes')
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('gc content/valid aligned bases')
#plt.legend()
plt.savefig(f'{output_folder}/c_gc_content.png', format='png',bbox_inches='tight')

plt.figure()
plt.title('Cancer Genes')
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_non_syn_percentage[speci], label=speci, s = size)
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('M_n/valid aligned bases')
#plt.legend()
plt.savefig(f'{output_folder}/c_M_n_content.svg', format='svg',bbox_inches='tight')


plt.figure()
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_GC_ending_percentage[speci], label=speci, s = size)
plt.title('Cancer Genes')
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('GC_ending/total 4 fold degen sites')
#plt.legend()
plt.savefig(f'{output_folder}/c_gc_ending_4fold.svg', format='svg',bbox_inches='tight')


