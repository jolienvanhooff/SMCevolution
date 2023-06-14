#!/usr/bin/env python3

import pandas as pd

euk5table = pd.read_csv('../phylogenetic_profiles.csv', index_col="Abbreviation")
species_abbrevs = euk5table.index.tolist()

with open('iTOL_label_names.txt', 'w') as f:
    f.write('LABELS\n')
    f.write('SEPARATOR TAB\n')
    f.write('DATA\n')
    for abbrev in species_abbrevs:
        species_name = euk5table.loc[abbrev]['Scientific name']
        f.write(abbrev + '\t' + species_name + '\n')
