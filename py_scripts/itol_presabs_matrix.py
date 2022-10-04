import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--complex', required=True)
parser.add_argument('-p', '--protein', nargs = '+', required=True)
parser.add_argument('-col', '--colour', required=True)

args = parser.parse_args()
args = dict(vars(args))
protein_complex = args.get('complex')
proteins = list(args.get('protein'))
colour = args.get('colour')

##Update pad naar de eukarya_species.tsv

eukarya_species = pd.read_csv('/home/max/NOBINFBACKUP/eukarya5_max/data_set/eukarya_species.tsv', sep = '\t', header = 'infer')
eukarya_species = eukarya_species['Abbreviation'].tolist()

field_labels = list()
for protein in proteins:
    field_labels.append(protein)
field_labels = '\t'.join(field_labels)

field_colours = '{}\t'.format(colour) * len(proteins)
print(len(proteins))
field_shapes = '1\t' * len(proteins)

##Update paden naar 1) waar je je output file wil hebben en 2) waar je lijstjes met ortholog ids staan

with open('itol_presabs/itol_presabs_{}.txt'.format(protein_complex), 'w') as f:
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
    for species in eukarya_species:
        species_dict = {}
        for protein in proteins:
            species_dict[protein] = '0'
            ids = pd.read_csv('all_ids/{}_all_ids.csv'.format(protein), sep = '\t', header = None)
            ids = ids[0].tolist()
            for i, seq_id in enumerate(ids):
                if '_' in seq_id:
                    seq_id_isolated = seq_id.split('_')[1]
                    ids[i] = seq_id_isolated
            ids = [x[0:6] for x in ids]
            if species in ids:
                species_dict[protein] = '1'
        protein_values = '\t'.join(species_dict.values())
        f.write(species + '\t' + protein_values + '\n')
