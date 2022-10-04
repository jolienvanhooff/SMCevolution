import pandas as pd

species_taxonomy = pd.read_csv('./eukarya_species.tsv', sep = '\t', header = 'infer')
species_list = species_taxonomy['Scientific name'].tolist()
species_taxonomy = species_taxonomy.set_index('Scientific name')

colour_codes = {"Stramenopiles" : "#750C19", "Alveolata" : "#762A83", "Rhizaria" : "#B85A07", "Telonemia" : "#74ADD1", "Haptista" : "#40004B", "Ancoracysta" : "#C2A5CF", "Chloroplastida" : "#1B7837", "Glaucophyta" : "#6BC987", "Rhodophyta" : "#F7001E", "Rhodelphis" : "#F26374", "Cryptista" : "#35978F", "Discoba" : "#C1DE3E", "Obazoa" : "#F59A22", "Amoebozoa" : "#F04695", "Ancyromonadida" : "#AB688A", "Metamonada" : "#3288BD", 'CRuMs' : "#1824AB", "Malawimonadidae" : "#5E5E5E",  "Hemimastigophora" : "#000000", "Picozoa" : "#2F6140", 'Breviatea' : '#f59a22', 'Apusomonadida' : '#F59A22', 'Opisthokonta' : '#F59A22'}

with open('iTOL_Taxonomy_BranchColour.txt', 'w') as f:
    f.write('DATASET_STYLE\n')
    f.write('SEPARATOR TAB\n')
    f.write('DATASET_LABEL\tTaxonomy\n')
    f.write('COLOR\t#000000\n')
    f.write('LEGEND_TITLE\tTaxonomy\n')
    f.write('LEGEND_SHAPES\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\n')
    f.write('LEGEND_COLORS\t#750C19\t#762A83\t#B85A07\t#74ADD1\t#40004B\t#C2A5CF\t#1B7837\t#6BC987\t#F7001E\t#F26374\t#35978F\t#C1DE3E\t#F59A22\t#F04695\t#AB688A\t#3288BD\t#1824AB\t#5E5E5E\t#000000\t#2F6140\t#f59a22\t#F59A22\t#F59A22\n')
    f.write('LEGEND_LABELS\tStramenopiles\tAlveolata\tRhizaria\tTelonemia\tHaptista\tAncoracysta\tChloroplastida\tGlaucophyta\tRhodophyta\tRhodelphis\tCryptista\tDiscoba\tObazoa\tAmoebozoa\tAncyromonadida\tMetamonada\tCRuMs\tMalawimonadidae\tHemimastigophora\tPicozoa\tBreviatea\tApusomonadida\tOpisthokonta\n')
    f.write('DATA\n')
    for species in species_list:
        taxonomy = species_taxonomy.loc[species]['relevant taxonomy']
        species = species.replace(' ', '_')
        colour_code = colour_codes.get(taxonomy)
        f.write('{}\tbranch\tclade\t{}\t1\tnormal\n'.format(species, colour_code))


