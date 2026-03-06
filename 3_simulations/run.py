from MutationTableGenerator import MutationTableGenerator
from CalculateSimulatedTarget import CalculateSimulatedTarget
from MultiplyTargets import MultiplyTargets
from MutationMatrixGenerator import MutationMatrixGenerator
from concurrent.futures import ProcessPoolExecutor
import sys
import pandas as pd
from NeutralEvolutionSimulator import NeutralEvolutionSimulator
from pathlib import Path
import shutil

#parameters to edit on each run
#target='hg38'
#anc='Anc4'
#N=200000

target='Anc4'
anc='Anc3'
N=150000

#target='fullTreeAnc212'
#anc='fullTreeAnc213'
#N=150000

#target='Anc1'
#anc='Anc0'
#N=170000

#target='Anc3'
#anc='Anc1'  #ancestor species
#N=300000 #effective population size
Nb=0.07  #modulates gc bias strength

#prarmeters keep constant across runs
episode = f'{target}_{anc}'
#input parameters
dir = f'/home/maria/run_simulations_cactus_clean/{episode}'
exon_file=f'{dir}/bed_files/exons_nodup_nocpg_flank.bed'  
gene_annotations='/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'  #gene annotations file #merged ids sorted 
mutation_prob_file=f'{dir}/intron_matrix/rates'  #mutation prob file from introns


nonsyn_data_prep=f'/home/maria/run_simulations_cactus_clean/scripts/2_run_cbase/{episode}/{episode}/{anc}->{target}/nonsyn_output_data_preparation_{episode}.txt'
syn_data_prep=f'/home/maria/run_simulations_cactus_clean/scripts/2_run_cbase/{episode}/{episode}/{anc}->{target}/syn_output_data_preparation_{episode}.txt'

possible_species = [target, anc]
auxiliary_dir = f'{dir}/sim_auxiliary'
target_dir = f'{dir}/sim_target_sizes'
output_dir= f'{dir}/output' #output dir


#already ran first part for episode 

print('generating mutation tables for ancestor species')
#create mutations table for ancestor species
generator = MutationTableGenerator(anc, possible_species, exon_file, gene_annotations, auxiliary_dir,mutation_prob_file,N, Nb)
generator.run()




#call anc target size
matrix_generator = MutationMatrixGenerator('', '', '', '', '')
print('calculating ancestor target')
target_caller = CalculateSimulatedTarget(
matrix_generator,
input_file = f'{auxiliary_dir}/{anc}_processed_gene_seqs.pkl',
output_dir = f'{target_dir}/{anc}')
target_caller.run_parallel()


def process_one_signature(sig) -> str:
    calc = MultiplyTargets(signature=sig,
                            input_glob= f'{target_dir}/{anc}/*',
                            output_dir=f"{output_dir}/{anc}",
                            gene_strand_file='/home/maria/filter_transcripts/output/exon_merged_ids_strands')
    calc.calc_for_sim()
    return 

DEFAULT_SIGNATURE_LIST = [
    "SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7a","SBS7b","SBS7c","SBS7d",
    "SBS8","SBS9","SBS10a","SBS10b","SBS10c","SBS10d","SBS11","SBS12","SBS13","SBS14",
    "SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19","SBS20","SBS21","SBS22a","SBS22b",
    "SBS23","SBS24","SBS25","SBS26","SBS28","SBS29","SBS30","SBS31","SBS32","SBS33",
    "SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40a","SBS40b","SBS40c","SBS41",
    "SBS42","SBS44","SBS84","SBS85","SBS86","SBS87","SBS88","SBS89","SBS90","SBS91",
    "SBS92","SBS93","SBS94","SBS96","SBS97","SBS98","SBS99"
]

print('mulitplying by signatures for ancestor species')
with ProcessPoolExecutor(max_workers=10) as executor:
    executor.map(process_one_signature, DEFAULT_SIGNATURE_LIST)


print('calculating mutation target for target species')
#call target size for target species
generator = MutationMatrixGenerator(target, possible_species, exon_file, gene_annotations, f'{target_dir}/{target}')
generator.run()

#call multiply by signatures for target species
def process_one_signature(sig) -> str:
    calc = MultiplyTargets(signature=sig,
                            input_glob= f'{target_dir}/{target}/*',
                            output_dir=f'{output_dir}/{target}',
                            gene_strand_file='/home/maria/filter_transcripts/output/exon_merged_ids_strands')
    calc.calc_for_sim()
    return 

DEFAULT_SIGNATURE_LIST = [
    "SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7a","SBS7b","SBS7c","SBS7d",
    "SBS8","SBS9","SBS10a","SBS10b","SBS10c","SBS10d","SBS11","SBS12","SBS13","SBS14",
    "SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19","SBS20","SBS21","SBS22a","SBS22b",
    "SBS23","SBS24","SBS25","SBS26","SBS28","SBS29","SBS30","SBS31","SBS32","SBS33",
    "SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40a","SBS40b","SBS40c","SBS41",
    "SBS42","SBS44","SBS84","SBS85","SBS86","SBS87","SBS88","SBS89","SBS90","SBS91",
    "SBS92","SBS93","SBS94","SBS96","SBS97","SBS98","SBS99"
]

print('mulitplying by signatures for target species')
with ProcessPoolExecutor(max_workers=10) as executor:
    executor.map(process_one_signature, DEFAULT_SIGNATURE_LIST)


print('processing cbase data')
# prcoess cbase data
syn_data_df = pd.read_csv(syn_data_prep, index_col=0, delimiter='\t')
nonsyn_data_df = pd.read_csv(nonsyn_data_prep, index_col=0, delimiter='\t')
data_df = pd.merge(left = syn_data_df, right=nonsyn_data_df, how= 'inner', on = ['gene', 'l_m',	'l_k', 'l_s',	'm_obs',	'k_obs',	's_obs',	'L_gene','N_samples=1'])
data_df['muts'] = data_df['s_obs'] + data_df['m_obs'] + data_df['k_obs']

cut_df = data_df[['lambda_n','lambda_s']]
cut_df.to_csv(f'{auxiliary_dir}/cbase.csv')

N_SUBS=data_df['muts'].sum()

#printing parameters
print(f'number of mutations {data_df['muts'].sum()}')
print(f'number of unique genes {len(cut_df.index.unique())}')
mut_df = data_df[data_df['muts']>0]
print(f'number of mutated genes {len(mut_df.index.unique())}')




#run neutral evolution simulator
sim_start_no = 0
no_of_sims = 100
def process_one_signature(args) -> str:
    sig,x = args
    calc = MultiplyTargets(signature=sig,
                            input_glob= f'{target_dir}/sim_{x}/*',
                            output_dir=f"{output_dir}/sim_{x}",
                            gene_strand_file='/home/maria/filter_transcripts/output/exon_merged_ids_strands')
    calc.calc_for_sim()
    return 



DEFAULT_SIGNATURE_LIST = [
    "SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7a","SBS7b","SBS7c","SBS7d",
    "SBS8","SBS9","SBS10a","SBS10b","SBS10c","SBS10d","SBS11","SBS12","SBS13","SBS14",
    "SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19","SBS20","SBS21","SBS22a","SBS22b",
    "SBS23","SBS24","SBS25","SBS26","SBS28","SBS29","SBS30","SBS31","SBS32","SBS33",
    "SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40a","SBS40b","SBS40c","SBS41",
    "SBS42","SBS44","SBS84","SBS85","SBS86","SBS87","SBS88","SBS89","SBS90","SBS91",
    "SBS92","SBS93","SBS94","SBS96","SBS97","SBS98","SBS99"
]
matrix_generator = MutationMatrixGenerator('', '', '', '', '')

for x in range(sim_start_no,no_of_sims):
    print(f'running simulation {x}')
    simulator = NeutralEvolutionSimulator(
        n_subs=N_SUBS,
        start_speci_file=f'{auxiliary_dir}/{anc}_processed_gene_seqs.pkl',
        cbase_output=f'{auxiliary_dir}/cbase.csv',
        syn_probs_dir=f'{auxiliary_dir}/syn_target',
        non_syn_probs_dir=f'{auxiliary_dir}/nonsyn_target',
        output=f'{output_dir}/simulated_genes_{x}.pkl',
        gene_name_map = '/home/maria/filter_transcripts/output/gene_name_id.csv')
    simulator.run()
    print(f'calculating targets for simulation {x}')
    target_caller = CalculateSimulatedTarget(
    matrix_generator,
    input_file = f'{output_dir}/simulated_genes_{x}.pkl',
    output_dir = f'{target_dir}/sim_{x}',
    )
    target_caller.run_parallel()
    print(f'mutiplying targets for simulation {x}')
    args_list = [(sig, x) for sig in DEFAULT_SIGNATURE_LIST]
    with ProcessPoolExecutor(max_workers=15) as executor:
        executor.map(process_one_signature, args_list)
    p=Path(f'{target_dir}/sim_{x}')
    if p.exists() and p.is_dir():
        shutil.rmtree(p)


