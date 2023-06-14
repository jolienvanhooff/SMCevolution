#!/usr/bin/env python3

import pandas as pd
import argparse
from smc_variables import complex_codes, complex_members, complex_colours

parser = argparse.ArgumentParser(description="make a dataset to upload to iTol for an SMC protein complex")
parser.add_argument('-c', '--complex', required=True, help="name of the protein complex: Cond, CondI, CondII, SMC56 or Coh")
parser.add_argument('-pe', '--proteinexclude', nargs = '+', required=False, help="name(s) of the protein(s) that are part of the complex but should be excluded from the dataset")
parser.add_argument('-pp', '--pprofile', default="../phylogenetic_profiles.csv", help="phylogenetics profile comma-separated table")

args = vars(parser.parse_args())
protein_complex_code = args.get('complex')
proteins_excluded = list(args.get('proteinexclude')) if args.get('proteinexclude') != None else []
eukarya_profiles = pd.read_csv(args.get('pprofile'), index_col='Abbreviation')
eukarya_species = eukarya_profiles.index.tolist()

# Get data from variables
protein_complex = complex_codes[protein_complex_code]
protein_complex_members = complex_members[protein_complex]
proteins = [p for p in protein_complex_members if p not in proteins_excluded]
complex_colour_fill = complex_colours[protein_complex][0]
complex_colour_stroke = complex_colours[protein_complex][1]

# Construct strings for writing formatting options to the output file
field_labels = list()
for protein in proteins:
    field_labels.append(protein)
field_labels = ','.join(field_labels)
field_colour_fill = ','.join([complex_colour_fill for n in range(0,len(proteins))])
field_shapes = ','.join(['1' for n in range(0,len(proteins))]) ##1:square 2:circle 3:star 4:right pointing triangle 5:left pointing triangle 6:checkmark   

# Generate dataset textfile (the output)
with open(f'itol_presabs_{protein_complex_code}.txt', 'w') as f:
    f.write('DATASET_BINARY\n')
    f.write('SEPARATOR COMMA\n')
    f.write('DATASET_LABEL,{}\n'.format(protein_complex))
    f.write('COLOR,{}\n'.format(complex_colour_fill))
    f.write('FIELD_LABELS,{}\n'.format(field_labels))
    f.write('FIELD_COLORS,{}\n'.format(field_colour_fill))
    f.write('FIELD_SHAPES,{}\n'.format(field_shapes))
    f.write('LEGEND_SHAPES,{}\n'.format(field_shapes))
    f.write('LEGEND_TITLE,{}\n'.format(protein_complex))
    f.write('LEGEND_COLORS,{}\n'.format(field_colour_fill))
    f.write('LEGEND_LABELS,{}\n'.format(field_labels))
    f.write('SHOW_LABELS,1\n')
    f.write('SIZE_FACTOR,0.8\n')
    f.write('LABEL_ROTATION,0\n')
    f.write('HEIGHT_FACTOR,1.1\n')
    f.write('SYMBOL_SPACING,0\n')
    #f.write('DASHED_LINES,1\n')
    f.write('MARGIN,-8\n')
    f.write('DATA\n')
    # Collect binary results (presence or absence) for each species
    for species in eukarya_species:
        protein_value = []
        for protein in proteins:
            protein_count = eukarya_profiles.loc[species, protein] 
            protein_value.append('1' if protein_count > 0 else '-1') # 0:empty, but lined shape in case a protein in absent in species -1: no symbol at all in absences
        protein_value_string = ','.join(protein_value)
        f.write(species + ',' + protein_value_string + '\n')
