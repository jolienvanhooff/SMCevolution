#!/usr/bin/env python3

import pandas as pd

species_df = pd.read_csv('../phylogenetic_profiles.csv', index_col="Abbreviation")
with open("../../euk5proteomes/euk5_tree_abbr.nwk", 'r') as t:
    tree = t.read()
species_tree_order = tree.replace("(", "").replace(")","").replace(";","").replace("\n", "")
species_ordered = species_tree_order.split(",")
species_df = species_df.reindex(species_ordered)

clades_colour_codes = {"Stramenopiles" : "#750C19", "Alveolata" : "#762A83", "Rhizaria" : "#B85A07", "Telonemia" : "#74ADD1", "Haptista" : "#40004B", "Ancoracysta" : "#C2A5CF", "Chloroplastida" : "#1B7837", "Glaucophyta" : "#6BC987", "Rhodophyta" : "#F7001E", "Rhodelphis" : "#F26374", "Cryptista" : "#35978F", "Discoba" : "#C1DE3E", "Obazoa" : "#F59A22", "Amoebozoa" : "#F04695", "Ancyromonadida" : "#AB688A", "Metamonada" : "#3288BD", 'CRuMs' : "#1824AB", "Malawimonadidae" : "#5E5E5E",  "Hemimastigophora" : "#000000", "Picozoa" : "#2F6140", 'Breviatea' : '#f59a22', 'Apusomonadida' : '#F59A22', 'Opisthokonta' : '#F59A22'}

## Dataset for clade colours (colours of branches)
with open('iTOL_clade_branchcolour.txt', 'w') as f:
    f.write('DATASET_STYLE\n')
    f.write('SEPARATOR COMMA\n')
    f.write('DATASET_LABEL,Clade\n')
    f.write('COLOR,#000000\n')
    f.write('LEGEND_TITLE,Clade\n')
    f.write('LEGEND_SHAPES,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1\n')
    f.write('LEGEND_COLORS,#750C19,#762A83,#B85A07,#74ADD1,#40004B,#C2A5CF,#1B7837,#6BC987,#F7001E,#F26374,#35978F,#C1DE3E,#F59A22,#F04695,#AB688A,#3288BD,#1824AB,#5E5E5E,#000000,#2F6140,#f59a22,#F59A22,#F59A22\n')
    f.write('LEGEND_LABELS,Stramenopiles,Alveolata,Rhizaria,Telonemia,Haptista,Ancoracysta,Chloroplastida,Glaucophyta,Rhodophyta,Rhodelphis,Cryptista,Discoba,Obazoa,Amoebozoa,Ancyromonadida,Metamonada,CRuMs,Malawimonadidae,Hemimastigophora,Picozoa,Breviatea,Apusomonadida,Opisthokonta\n')
    f.write('LEGEND_SHAPE_SCALES,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1\n')
    f.write('DATA\n')
    for clade, colour in clades_colour_codes.items():
        members = species_df[species_df['relevant taxonomy'] == clade].index.tolist()
        if len(members) > 1:
            branch = f"{members[0]}|{members[-1]}"
        elif len(members) == 1:
            branch = members[0]
        else:
            continue
        f.write(f"{branch},branch,clade,{colour},1,normal\n")

## Dataset for clade names
with open('iTOL_clade_label.txt', 'w') as l:
    l.write("LABELS\n")
    l.write("SEPARATOR COMMA\n")
    l.write("DATA\n")
    for clade in clades_colour_codes.keys():
        members = species_df[species_df['relevant taxonomy'] == clade].index.tolist()
        if len(members) > 1:
            branch = f"{members[0]}|{members[-1]}"
        elif len(members) == 1:
            branch = members[0]
        else:
            continue
        l.write(f"{branch},{clade},Supergroup\n") 
