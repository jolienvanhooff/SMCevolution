#!/usr/bin/env python3

import pandas as pd

euk5table = pd.read_csv('eukarya_species.tsv', sep = '\t', header = 'infer')
species_abbrevs = euk5table['Abbreviation'].tolist()
euk5table = euk5table.set_index('Abbreviation')

with open('iTOL_label_names.txt', 'w') as f:
    f.write('LABELS\n')
    f.write('SEPARATOR TAB\n')
    f.write('DATA\n')
    for abbrev in species_abbrevs:
        species_name = euk5table.loc[abbrev]['Scientific name']
        f.write(abbrev + '\t' + species_name + '\n')
