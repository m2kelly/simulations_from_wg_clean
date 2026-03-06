target='Anc4'
anc='Anc3'

episode="${target}_${anc}"
auxiliary="/home/maria/run_simulations_cactus_clean/${episode}/bed_files"
input="${auxiliary}/${target}_${anc}.bed"
cpg_coords="${auxiliary}/combined_cpgs.bed" #species specific cpg coords
aligned_file="${auxiliary}/exons_nodup_nocpg_flank_aligned.bed"
output_file="${auxiliary}/exons_nodup_nocpg_flank.bed"


dups='/home/maria/data/segmental_duplicates_cleaned.bed'
exon_coords='/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
exon_nooverlap='/home/maria/filter_transcripts/output/exon_filt.bed'
coords="${auxiliary}/exons_nodup_nocpg_flank_coords.bed" #intermediate file


tail -n +2 $input | \
cut -f1-3 | \
bedtools intersect -a - -b $exon_nooverlap | \
cut -f1-3 | \
bedtools subtract -a - -b $dups | \
bedtools subtract -a - -b $cpg_coords \
> $coords

tail -n +2 $input | \
bedtools intersect -wa -wb -a $coords -b - | \
cut -f 1-3,5-8 | \
bedtools intersect -wa -wb -a - -b $exon_coords | \
uniq \
> $aligned_file
#[ 'chr','extract start', 'extract end', 'seq_start','seq_end', 'target_seq','seq', 'chr_copy','exon start', 'exon end', 'gene','gene_name','strand']

rm $coords
python3 /home/maria/run_simulations_cactus_clean/scripts/1_extract_exons/extract_exon_coords.py $aligned_file $output_file
rm $aligned_file

