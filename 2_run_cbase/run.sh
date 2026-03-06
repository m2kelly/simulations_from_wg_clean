#prepare file
target='Anc4'
anc='Anc3'

#ALSO NEED TO SET TARGET AND ANC IN INPUT.PY
episode="${target}_${anc}"

dir="/home/maria/run_simulations_cactus_clean/${episode}"
cpg_coords="${dir}/bed_files/combined_cpgs.bed" #species specific cpg coords

mkdir -p "$dir/cbase"

cat /home/maria/data/segmental_duplicates_cleaned.bed $cpg_coords \
  | bedtools sort -i - \
  | bedtools merge -i - \
  > $dir/cbase/cpg_dup.bed

python3 /home/maria/run_simulations_cactus_clean/scripts/2_run_cbase/edit_inputbed_to_gtf_cpg.py $dir/cbase/cpg_dup.bed $dir/cbase/cpg_dup.gtf 

rm $dir/cbase/cpg_dup.bed

python3 /home/maria/run_simulations_cactus_clean/scripts/2_run_cbase/main.py

