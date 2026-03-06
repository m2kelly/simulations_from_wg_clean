
from IntronMatrixNoCpG import IntronMatrixNoCpG
import sys

target=sys.argv[1]  #target species
anc=sys.argv[2]  #ancestor species
intron_file=sys.argv[3]  #exon file
#test
'''
target='fullTreeAnc212'
anc='fullTreeAnc213'
intron_file='/home/maria/run_simulations_cactus_clean/fullTreeAnc212_fullTreeAnc213/bed_files/introns_nodup_nocpg_flank_aligned.bed'
'''
episode = f'{target}_{anc}'
output_dir = f'/home/maria/run_simulations_cactus_clean/{episode}/intron_matrix'


species = ''
target_species = ''
possible_species = []


introns=IntronMatrixNoCpG(intron_file, output_dir,species, target_species, possible_species)
introns.run()
