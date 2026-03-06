import pandas as pd
import sys


in_file  = sys.argv[1]
out_file = sys.argv[2]

#in_file  = '/home/maria/filter_transcripts/output/exon_merged_ids.bed'


# Input has: chr, start, end, gene_id, gene, strand
#switching gene and gene_id, to use gene_id as gene name 
df = pd.read_csv(in_file, sep='\t',
                 names=['chr', 'start', 'end'],
                 dtype={'chr': str, 'start': int, 'end': int})

# If your input is BED-like (0-based, half-open) and you want GTF (1-based, closed), uncomment this:
df['start'] = df['start'] + 1

# Required GTF/GFF columns
df['source']  = 'refGene'          # your chosen source
df['feature'] = 'merged_prohibited'    # your chosen feature/type
df['score']   = '.'                # no score
df['frame']   = '.'                # no frame
df['strand'] = '+' #should not effect?
# GTF-style attributes (you can add more key-value pairs if needed)
#df['attributes'] = 'gene_id "' + df['gene_id'] + '"; gene_name "' + df['gene'] + '";'
#df['attributes'] = 'gene_name "' + df['gene'] + '";'
df['attributes'] = ''

# Reorder to 9 columns (GTF order)
df_out = df[['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame','attributes']]

# Save
df_out.to_csv(out_file, sep='\t', header=False, index=False)
print(f"Written: {out_file}")