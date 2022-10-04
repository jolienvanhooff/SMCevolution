import pandas as pd

busco_output = pd.read_csv('BUSCO_output.tsv', header = 'infer', sep = '\t')
species_list = busco_output['Species'].tolist()
busco_output = busco_output.set_index('Species')

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
        f.write(species + '\t-1\t20\t' + str(busco_output.loc[species]['Percentage.BUSCO.complete']) + '\t' + str(busco_output.loc[species]['Percentage.fragmented']) + '\t' + str(busco_output.loc[species]['Percentage.missing']) + '\n')
