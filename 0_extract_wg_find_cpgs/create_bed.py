'''
input- cactus hal2maf output
output - bed file of human coords and coorepsonding alligned seqs
pipeline
1) parse all lines, species lines begins with s, and every alignment bloc starts with a
2) only parse and save lines with species in species_output
2) at end of each alignememnt bloclk if exactly one line per species then save line in bed format with human coordinates
3) for each row in bed file remove any insertions in human seq, and remove correposning bases in other seqs 

if do not want to output ref species, but want to use for alignement and removing insertions
then add tag include ref species = False (but still put it in species output)
'''


import gzip
import pandas as pd
import sys
import os 

target=sys.argv[1]
anc=sys.argv[2]
input = sys.argv[3]
output_dir = sys.argv[4]


ref_species = 'hg38'
#ref_species = 'Homo_sapiens'

species_output = [target, anc]
if ref_species in species_output:
    species_extract = species_output
else:
    species_extract = species_output + [ref_species]
print(species_extract)
output = f'{output_dir}/{target}_{anc}.bed'
#make outputdir
mkdir_cmd = f'mkdir -p {output_dir}/'
os.system(mkdir_cmd)


sep = "\t"

VALID_CHR = {f'chr{x}' for x in range(1, 23)}
BASES = set("ACGT")

def parse_line(line):
    fields = line.rstrip('\n').split(sep)
    species = fields[1].split(".")[0]
    if species == ref_species:
        chr = fields[1].split(".")[1]
        strand = fields[4]
        start_pos = int(fields[2])
        length = int(fields[3])
        end_pos = start_pos + length
        sequence = fields[6]
        return species, [sequence, length, chr, start_pos, end_pos, strand]
    else:
        sequence = fields[6]
        length = int(fields[3])
        return species,  [sequence, length]

def save_line(species, parsed_lines, coords):
    line_dict = {}
    line_dict['chr'] = coords[0]
    line_dict['start'] = coords[1]
    line_dict['end'] = coords[2]
    for i, speci in enumerate(species):
        seq = parsed_lines[i][0]
        line_dict[speci] = seq
    return line_dict


#editting to use one boolean mask
def remove_insertions_in_human(row):
    hg = row[ref_species]
    if '-' in hg:
        keep = [c != '-' for c in hg]   # boolean mask once
        filt_speci={}
        for speci in species_extract:
            no_hum_ins = ''.join(b for k, b in zip(keep, row[speci]) if k)
            filt_speci[speci] = ''.join('N' if b =='-' else b for b in no_hum_ins) #mask with N so works with cpg masker
    else:
        filt_speci = {speci: ''.join('N' if b =='-' else b for b in row[speci]) for speci in species_extract}

    # zip(mask, seq) is fast; avoids membership checks
    return pd.Series(filt_speci)


# -------------------------
# NEW: batch flushing helper
# -------------------------
BATCH_SIZE = 1000  # adjust up/down depending on memory/speed
first = True

def flush_batch(all_lines):
    """Convert current batch to df, run same transforms, append to output, then clear."""
    global first
    if not all_lines:
        return

    exon_df = pd.DataFrame(all_lines)
    
    # processing steps
    exon_df[species_extract] = exon_df.apply(lambda row: remove_insertions_in_human(row), axis=1)
    exon_df = exon_df[['chr','start','end']+species_output]
    

    # append to disk, ensure col order always the same 
    exon_df[['chr','start','end']+species_output].to_csv(
        output,
        index=False,
        sep='\t',
        mode='a',
        header=(first)
    )
    first=False
    all_lines.clear()
    return

# -------------------------
# streaming parse (same logic, just flush periodically)
# -------------------------
all_lines = []
print('parsing lines in maf file')

with gzip.open(input, 'rt') as f_in:
    next(f_in); next(f_in); next(f_in)

    species = []
    parsed_lines = []
    coords = []

    for line in f_in:
    #TESTING NF
    #for line in islice(f_in,100000000):
        if line[0] == 'a':
            
            uniq_species = list(set(species))
            if len(species) == len(uniq_species):  # no repeated seqs
                if set(species) == set(species_extract):
                    #save lines if valid chr
                    if coords[0] in VALID_CHR:
                        line_dict = save_line(species, parsed_lines, coords)
                        all_lines.append(line_dict)

                    # flush every BATCH_SIZE valid blocks
                    if len(all_lines) >= BATCH_SIZE:
                        #print(f'flushing batch of {len(all_lines)}')
                        flush_batch(all_lines)
                        #break
                        

            # reset for next block
            parsed_lines = []
            species = []
            coords = []

        elif line[0] == 's':

            speci, parsed_line = parse_line(line)
            
            if speci == ref_species:
                coords = parsed_line[2:]
                parsed_lines.append(parsed_line[:2])
                species.append(speci)
            elif speci in species_output:
                parsed_lines.append(parsed_line)
                species.append(speci)

    # save last block if valid
    uniq_species = list(set(species))
    if len(species) == len(uniq_species): #if no duplicates 
        if set(species) == set(species_extract):
            line_dict = save_line(species, parsed_lines, coords)
            all_lines.append(line_dict)

# NEW: flush remainder
print(f'flushing final batch of {len(all_lines)}')
flush_batch(all_lines)

print("done")
