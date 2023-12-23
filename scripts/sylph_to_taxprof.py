import pandas as pd
from collections import defaultdict
import sys
import csv
import gzip
import argparse

def genome_file_to_gtdb_acc(file_name):
    return(file_name.split('/')[-1].split("_genomic")[0])

parser = argparse.ArgumentParser(description='Convert a sylph output tsv file to a set of output files in MetaPhlAn format.')
parser.add_argument("-s", "--sylph", help="sylph output file", type=str, required=True)
parser.add_argument("-o", "--output-prefix", help="prefix of the outputs. Output files will be prefix + Sample_file +.sylphmpa", type=str, default = "")
parser.add_argument("-m", "--metadata", help = "metadata file converting genome files to taxnomic identifiers.", type = str, required=True) 
args = parser.parse_args()

# Read sylph's output TSV file into a Pandas DataFrame
df = pd.read_csv(args.sylph, sep='\t')

### This is a dictionary that contains the genome_file 
### to taxonomy string mapping. It should be like
### {'my_genome.fna.gz' : b__Bacteria;...}
genome_to_taxonomy = dict()

### Process gzip file instead if extension detected
if '.gz' in args.metadata or '.gzip' in args.metadata:
    f=gzip.open(args.metadata,'rt')
else:
    f = open(args.metadata,'r')

### Tag each taxonomy string with a t__ strain level identifier
for row in f:
    spl = row.rstrip().split();
    accession = spl[0]
    taxonomy = ' '.join(spl[1:]).rstrip() + ';t__' + accession
    genome_to_taxonomy[accession] = taxonomy

### Group by sample file. Output one file fo reach sample file. 
grouped = df.groupby('Sample_file')

for sample_file, group_df in grouped:
    out = sample_file.split('/')[-1]
    out_file = args.output_prefix + out + '.sylphmpa'
    of = open(out_file,'w')

    tax_abundance = defaultdict(float)
    seq_abundance = defaultdict(float)
    ani_dict = defaultdict(float)

    # Iterate over each row in the DataFrame
    for idx, row in group_df.iterrows():

        # Parse the genome file... assume the file is in gtdb format.
        # This can be changed. 
        genome_file = genome_file_to_gtdb_acc(row['Genome_file'])
        ani = float(row['Adjusted_ANI'])

        if genome_file in genome_to_taxonomy:
            tax_str = genome_to_taxonomy[genome_file]
        elif genome_file +'.gz' in genome_to_taxonomy:
            tax_str = genome_to_taxonomy[genome_file + '.gz']
        else:
            tax_str = 'NO_TAXONOMY;t__' + genome_file

        abundance = float(row['Sequence_abundance'])
        rel_abundance = float(row['Taxonomic_abundance'])

        # Split taxonomy string into levels and update abundance
        tax_levels = tax_str.split(';')
        cur_tax = ''
        for level in tax_levels:
            if cur_tax:
                cur_tax += '|'
            cur_tax += level
            tax_abundance[cur_tax] += rel_abundance
            seq_abundance[cur_tax] += abundance
            if 't__' in cur_tax:
                ani_dict[cur_tax] = ani

    # Print the CAMI BioBoxes profiling format
    of.write(f"#SampleID\t{sample_file}\n")
    of.write("#clade_name\trelative_abundance\tsequence_abundance\tANI (if strain-level)\n")

    level_to_key = dict()
    for key in tax_abundance.keys():
        level = len(key.split('|'))
        if level not in level_to_key:
            level_to_key[level] = [key]
        else:
            level_to_key[level].append(key)

    sorted_keys = sorted(level_to_key.keys())

    for level in sorted_keys:
        keys_for_level = sorted(level_to_key[level], key = lambda x: tax_abundance[x], reverse=True)
        for tax in keys_for_level:
            if tax in ani_dict:
                of.write(f"{tax}\t{tax_abundance[tax]}\t{seq_abundance[tax]}\t{ani_dict[tax]}\n")
            else:
                of.write(f"{tax}\t{tax_abundance[tax]}\t{seq_abundance[tax]}\tNA\n")

