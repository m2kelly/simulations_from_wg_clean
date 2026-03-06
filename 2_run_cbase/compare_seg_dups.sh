
#prepare file
cat /home/maria/data/segmental_duplicates_cleaned.bed /home/maria/data/cpgIslandExt.bed \
  | bedtools sort -i - \
  | bedtools merge -i - \
  > /home/maria/hossam_cbase_exons/cpg_dup.bed


#checking things

#count total percent of bases removed by cpgs and dups 
cut -f1-3 /home/maria/filter_transcripts/output/exon_merged_ids.bed | \
bedtools intersect -a - -b /home/maria/hossam_cbase_exons/cpg_dup.bed | \
uniq | \
# BED (tab-separated, 0-based, half-open)
awk 'BEGIN{OFS="\t"} {len=$3-$2; if(len>0) sum+=len} END{print sum+0}' 
#5704079
#removing 17% of bases by removing cpg and dups

cut -f1-3 /home/maria/filter_transcripts/output/exon_merged_ids.bed \ | 
awk 'BEGIN{OFS="\t"} {len=$3-$2; if(len>0) sum+=len} END{print sum+0}'
#33634986

cut -f1-3 /home/maria/filter_transcripts/output/exon_merged_ids.bed \ | 
bedtools intersect -a - -b /home/maria/data/segmental_duplicates_cleaned.bed | \
uniq | \
# BED (tab-separated, 0-based, half-open)
awk 'BEGIN{OFS="\t"} {len=$3-$2; if(len>0) sum+=len} END{print sum+0}' 
#2848071
#removing 8.5 percent from segmental duplicates 

cat /home/maria/data/segmental_duplicates_cleaned.bed > test
cat /home/maria/data/cpgIslandExt.bed >> test

cut -f1,4-5 /home/maria/hossam_generate_intron/Input/seg_dubsandcpg_islands.gtf | \
diff - test | \
grep 'chr1' | \
head -n 100


cut -f1,4-5 /home/maria/hossam_generate_intron/Input/seg_dubsandcpg_islands.gtf > test
bedtools intersect -wa -wb -a test -b test | \
awk -F'\t' 'BEGIN{OFS="\t"} ($2!=$5) {print}'


bedtools intersect -wa -wb -a /home/maria/hossam_cbase_exons/cpg_dup.bed -b /home/maria/hossam_cbase_exons/cpg_dup.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} ($2!=$5) {print}'


