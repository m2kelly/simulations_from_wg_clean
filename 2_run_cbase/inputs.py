#inputs chnage per run
target='Anc4'
anc='Anc3'
#maf_file='/home/maria/effect_pop_bins/data/wg_otter.maf.gz'

maf_file  ='/home/maria/effect_pop_bins/data/homo_8_primate.maf.gz'
coord_spcs = 'hg38'
#coord_spcs = 'Homo_sapiens'
episode = f'{target}_{anc}'
excluded_regions = {f"/home/maria/run_simulations_cactus_clean/{episode}/cbase/cpg_dup.gtf": ['merged_prohibited']}

#need to chnage reference as hg38 or Homo_sapiens
external_matrix = f'/home/maria/run_simulations_cactus_clean/{episode}/intron_matrix/rates'

output_name = episode
output_directory = episode

coding = False
#analysis_pairs = [{'target':'hg38', 'reference':'Anc4'},{'target':'Anc4', 'reference':'Anc3'}]
analysis_pairs = [{'target':target, 'reference':anc}]




annotations_file = '/home/maria/hossam_cbase_exons/exon_merged_ids.gtf'
#must be 1 based gtf file 

region_coding_dict = {'first_coding_exon':True, 'internal_exon':True,
                      'last_coding_exon':True, 'one_exon':True,
                      'first_of_twos':True, 'second_of_twos':True,}



groups_to_collapse = [ ['first_coding_exon', 'first_of_twos',
                       'last_coding_exon', 'second_of_twos',
                       'internal_exon', 'one_exon']]


collapse_threshold = 999999999 #minimim number of mutations in the matrix, otherwise collapse

gene_min_size = 100




'''analysis_pairs = [{'reference':'HCLCA', 'target':'homo_sapiens'},
                  {'reference':'HCGLCA', 'target':'HCLCA'},
                  {'reference':'HCGOrLCA', 'target':'HCGLCA'},
                  {'reference':'HCGOrGib', 'target':'HCGOrLCA'},
                  {'reference':'HCGLCA', 'target':'homo_sapiens'},
                  {'reference':'HCGOrGib', 'target':'HCGLCA'},
                  ]'''

