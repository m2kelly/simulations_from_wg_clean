from MultiplyTargets import MultiplyTargets
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

#parameters to edit on each run

anc='Anc3'
target='Anc4'


#prarmeters keep constant across runs
episode = f'{target}_{anc}'
#input parameters
dir = f'/home/maria/run_simulations_cactus_clean/{episode}'

possible_species = [target, anc]
auxiliary_dir = f'{dir}/sim_auxiliary'
target_dir = f'{dir}/sim_target_sizes'
output_dir= f'{dir}/output_intron_matrix' #output dir
mut_matrix_file=f'{dir}/intron_matrix/rates'

mut_matrix=pd.read_csv(mut_matrix_file,index_col=0,sep='\t')
mut_matrix=mut_matrix/mut_matrix.sum().sum()

def process_one_signature(sig) -> str:
    calc = MultiplyTargets(signature=sig,
                            input_glob= f'{target_dir}/{anc}/*',
                            output_dir=f"{output_dir}/{anc}",
                            gene_strand_file='/home/maria/filter_transcripts/output/exon_merged_ids_strands')
    calc.sig='intron'
    calc.sig_df_neg=mut_matrix

    calc.sig_df_pos=mut_matrix
    calc.calc_for_sim()

    return 

DEFAULT_SIGNATURE_LIST = [
    "SBS5"
]

print('mulitplying by signatures for ancestor species')
with ProcessPoolExecutor(max_workers=10) as executor:
    executor.map(process_one_signature, DEFAULT_SIGNATURE_LIST)



def process_one_signature(sig) -> str:
    calc = MultiplyTargets(signature=sig,
                            input_glob= f'{target_dir}/{target}/*',
                            output_dir=f"{output_dir}/{target}",
                            gene_strand_file='/home/maria/filter_transcripts/output/exon_merged_ids_strands')
    calc.sig='intron'
    calc.sig_df_neg=mut_matrix

    calc.sig_df_pos=mut_matrix
    calc.calc_for_sim()

    return 

DEFAULT_SIGNATURE_LIST = [
    "SBS5"
]

print('mulitplying by signatures for target species')
with ProcessPoolExecutor(max_workers=10) as executor:
    executor.map(process_one_signature, DEFAULT_SIGNATURE_LIST)

