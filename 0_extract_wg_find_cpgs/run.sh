#not filtering for mislaigned in both species here-should be done for calculating target sizes
#this script creates bed files with cpg coords for both species in alignment
#then afterwards do th efiltering step for only alligned regions in both species when getting 
#intron and exon coords 
target='Anc4'
anc='Anc3'
maf='/home/maria/effect_pop_bins/data/homo_8_primate.maf.gz'

#maf='/home/maria/effect_pop_bins/data/wg_otter.maf.gz'


episode="${target}_${anc}"
auxiliary="/home/maria/run_simulations_cactus_clean/${episode}/bed_files"
bed="${auxiliary}/${target}_${anc}.bed"
cpg_coords=$auxiliary/combined_cpgs.bed

#make directory
mkdir -p "$auxiliary"

#need to run before removing letters that do not align in both 
#so make bed and remove insertions
python3 /home/maria/run_simulations_cactus_clean/scripts/0_extract_wg_find_cpgs/create_bed.py $target $anc $maf $auxiliary

#then create fasta file and feed to cpg_lh 
#this outputs bedfile with cpg coords
#run for both species


tail -n+2 $bed | cut -f1-3,4  | \
awk 'BEGIN{OFS=""} {print ">"$1":"$2"-"$3"\n"$4}' | ./cpg_lh /dev/stdin | \
awk '{$2 = $2 - 1; width = $3 - $2;
   printf("%s\t%d\t%s\t%s %s\t%s\t%s\t%0.0f\t%0.1f\t%s\t%s\n",
    $1, $2, $3, $5, $6, width, $6, width*$7*0.01, 100.0*2*$6/width, $7, $9);}'| \
sort -k1,1 -k2,2n  | \
awk 'BEGIN{OFS="\t"}{
  split($1,a,":");
  split(a[2],b,"-");

  chrom = a[1];
  winStart = b[1];

  islStart = $2;
  islEnd   = $3;

  bedStart = winStart + islStart - 2;
  bedEnd   = winStart + islEnd   - 1;

  print chrom, bedStart, bedEnd;
}' | sort -k1,1V -k2,2n > $auxiliary/target_cpgs.bed


tail -n+2 $bed | cut -f1-3,5  | \
awk 'BEGIN{OFS=""} {print ">"$1":"$2"-"$3"\n"$4}' | ./cpg_lh /dev/stdin | \
awk '{$2 = $2 - 1; width = $3 - $2;
   printf("%s\t%d\t%s\t%s %s\t%s\t%s\t%0.0f\t%0.1f\t%s\t%s\n",
    $1, $2, $3, $5, $6, width, $6, width*$7*0.01, 100.0*2*$6/width, $7, $9);}'| \
sort -k1,1 -k2,2n  | \
awk 'BEGIN{OFS="\t"}{
  split($1,a,":");
  split(a[2],b,"-");

  chrom = a[1];
  winStart = b[1];

  islStart = $2;
  islEnd   = $3;

  bedStart = winStart + islStart - 2;
  bedEnd   = winStart + islEnd   - 1;

  name = "CpG";
  score = $5;
  strand = ".";

  print chrom, bedStart, bedEnd;
}' | sort -k1,1V -k2,2n > $auxiliary/anc_cpgs.bed

#merge cpg bed files
cat $auxiliary/anc_cpgs.bed $auxiliary/target_cpgs.bed | sort -k1,1V -k2,2n | bedtools merge > $cpg_coords


#checl length of cpg regions
#awk '{sum += $3 - $2} END {print sum}' anc_cpgs.bed #22684944
#awk '{sum += $3 - $2} END {print sum}' target_cpgs.bed #22794372
#awk '{sum += $3 - $2} END {print sum}' combined_cpgs.bed #23695565
#awk '{sum += $3 - $2} END {print sum}' $bed #2686943817
#outer join of cpg islands-remove if island in either species
