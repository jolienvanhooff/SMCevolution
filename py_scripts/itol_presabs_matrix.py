#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="make a dataset to upload to iTol for a given protein complex")
parser.add_argument('-c', '--complex', required=True, help="name of the protein complex, e.g. Condensin I")
parser.add_argument('-p', '--protein', nargs = '+', required=True, help="name of the proteins constituting the complex and for which there is a phylogenetic profile in the phylogenetic profiles table")
parser.add_argument('-col', '--colour', required=True, type=str, help="colour for boxes indicating a protein's presence in a taxon; use quotation marks, e.g. '#FFC300'")
parser.add_argument('-pp', '--pprofile', default="../20220930_condensin_cohesin_SMC56.csv", help="phylogenetics profile comma-separated table")

args = parser.parse_args()
args = dict(vars(args))
protein_complex = args.get('complex')
proteins = list(args.get('protein'))
colour = args.get('colour')
eukarya_profiles = pd.read_csv(args.get('pprofile'), index_col='Abbreviation')
eukarya_species = eukarya_profiles.index.tolist()

# Construct strings for writing formatting options to the output file
field_labels = list()
for protein in proteins:
    field_labels.append(protein)
field_labels = '\t'.join(field_labels)
field_colours = '{}\t'.format(colour) * len(proteins)
field_shapes = '1\t' * len(proteins)

# Generate dataset textfile (the output)
with open(f'itol_presabs_{protein_complex}.txt', 'w') as f:
    f.write('DATASET_BINARY\n')
    f.write('SEPARATOR TAB\n')
    f.write('DATASET_LABEL\t{}\n'.format(protein_complex))
    f.write('COLOR\t{}\n'.format(colour))
    f.write('FIELD_LABELS\t{}\n'.format(field_labels))
    f.write('FIELD_COLORS\t{}\n'.format(field_colours))
    f.write('FIELD_SHAPES\t{}\n'.format(field_shapes))
    f.write('LEGEND_SHAPES\t{}\n'.format(field_shapes))
    f.write('LEGEND_TITLE\t{}\n'.format(protein_complex))
    f.write('LEGEND_COLORS\t{}\n'.format(field_colours))
    f.write('LEGEND_LABELS\t{}\n'.format(field_labels))
    f.write('DATA\n')
    # Collect binary results (presence or absence) for each species
    for species in eukarya_species:
        protein_value = []
        for protein in proteins:
            protein_count = eukarya_profiles.loc[species, protein] 
            protein_value.append('1' if protein_count > 0 else '0')
        protein_value_string = '\t'.join(protein_value)
        f.write(species + '\t' + protein_value_string + '\n')
