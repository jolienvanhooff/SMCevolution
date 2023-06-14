#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="make a dataset to upload to iTol to visualize BUSCO completeness scores")
parser.add_argument('-pp', metavar='pprofile', default="../phylogenetic_profiles.csv", help="phylogenetics profile comma-separated table")

args = parser.parse_args()
busco_output = pd.read_csv(args.pp, index_col='Abbreviation')
species_list = busco_output.index.tolist()

with open('iTOL_BUSCO_piechart.txt', 'w') as f:
    f.write('DATASET_PIECHART \n')
    f.write('SEPARATOR TAB\n')
    f.write('DATASET_LABEL\tBUSCO_score\n')
    f.write('COLOR\t#000000\n')
    f.write('FIELD_COLORS\t#000000\t#999999\t#ffffff\n')
    f.write('FIELD_LABELS\tComplete\tFragmented\tMissing\n')
    f.write('LEGEND_TITLE\tBUSCO_score\n')
    f.write('LEGEND_SHAPES\t2\t2\t2\n')
    f.write('LEGEND_COLORS\t#000000\t#999999\t#ffffff\n')
    f.write('LEGEND_LABELS\tComplete\tFragmented\tMissing\n')
    f.write('DATA\n')
    for species in species_list:
        f.write(species + '\t-1\t20\t' + str(busco_output.loc[species]['BUSCO_completeness']) + '\t' + str(busco_output.loc[species]['BUSCO_fragmented']) + '\t' + str(busco_output.loc[species]['BUSCO_missing']) + '\n')
