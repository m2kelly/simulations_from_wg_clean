
target='Anc4'
anc='Anc3'

episode="${target}_${anc}"
auxiliary="/home/maria/run_simulations_cactus_clean/${episode}/bed_files"
input="${auxiliary}/${target}_${anc}.bed"
cpg_coords="${auxiliary}/combined_cpgs.bed" #species specific cpg coords
aligned_file="${auxiliary}/introns_nodup_nocpg_flank_aligned.bed"

coords="${auxiliary}/introns_nodup_nocpg_flank_coords.bed" #intermediate file


bedtools subtract -a '/home/maria/filter_transcripts/output/intron_cords.bed' -b /home/maria/data/segmental_duplicates_cleaned.bed | \
bedtools subtract -a - -b $cpg_coords \
> $coords

tail -n +2 $input | \
bedtools intersect -wa -wb -a - -b $coords > $aligned_file 
rm $coords

python3 /home/maria/run_simulations_cactus_clean/scripts/1_extract_intron_matrix/main.py $target $anc $aligned_file
