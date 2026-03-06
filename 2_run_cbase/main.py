import glob
import sys
from CBaSE.interpreter import Interpreter
from CBaSE.models_fitter import ModelsFitter
from CBaSE.models_fitter_nonsyn import ModelsFitterNonsyn
from joblib import Parallel, delayed
import CBaSE.qvalues as qv
import inputs


#not fitting lambda_n-if want to do this use cbase_lambda_n - or simply chnage svae output_data_prep to save syn_outut_data_prep in
# file models_fitter line 170 
def analyze_interpreted(pair, inputs):
        source_seq =pair['reference'] ;target_seq = pair['target']
        ep_name = f'{source_seq}->{target_seq}'
        working_dir = f'{inputs.output_directory}/{inputs.output_name}/{ep_name}/'
        mks_file = f'{working_dir}/output_data_preparation_{inputs.output_name}.txt'
        ModelsFitterNonsyn(mks_file, working_dir, inputs.output_name)
        ModelsFitter(mks_file, working_dir, inputs.output_name)
        qv.compute_q_values(working_dir, inputs.output_name)
        
        exit()

if __name__ == "__main__":

    try: inputs.output_directory
    except: inputs.output_directory = 'Output'

    Interpreter(inputs)
    #jobs_no = min(8, len(inputs.analysis_pairs))
    Parallel(n_jobs=1, verbose=1)(map(delayed(analyze_interpreted), inputs.analysis_pairs,
                                                                         [inputs]*len(inputs.analysis_pairs)))

#hossam original cbase structure
'''
def analyze_interpreted(pair, inputs):
        source_seq =pair['reference'] ;target_seq = pair['target']
        ep_name = f'{source_seq}->{target_seq}'
        working_dir = f'{inputs.output_directory}/{inputs.output_name}/{ep_name}/'
        mks_file = f'{working_dir}/output_data_preparation_{inputs.output_name}.txt'
        ModelsFitter(mks_file, working_dir, inputs.output_name)
        qv.compute_q_values(working_dir, inputs.output_name)
        exit()
'''