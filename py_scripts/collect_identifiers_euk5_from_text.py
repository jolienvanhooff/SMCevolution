#!/usr/bin/env python3

import pandas as pd
import argparse
import os
from sys import argv
import re

parser = argparse.ArgumentParser(description="Gets the phyletic profiles of proteins from eukarya.v5 ('euk5') based on flat text (.txt) files with identifiers")
parser.add_argument("-t", metavar="table", type=str, default="../euk5proteomes/Euk5FinalSet.adjust.busco.euk5_tree_abbrev.csv", help="input table")
parser.add_argument("-d", metavar="base_directory", type=str, default="../proteins/", help="directory containing all protein folders")
parser.add_argument("-p", metavar="protein_order", type=str, default="SMC2,SMC4,CAPH,CAPG,CAPD2,CAPH2,CAPG2,CAPD3,SMC1,SMC3,Scc1,Rec8,Scc3,PDS5,NIPBL,MAU2,WAPL,Eco1,Securin,Sororin,Haspin,Shugoshin,Separase,CTCF,SMC5,SMC6,Nse4,Nse1,Nse3,Nse2,Nse5,Nse6", help="protein set to search for; ordered")
parser.add_argument("-o", metavar="output", type=str, default="phylogenetic_profiles_identifiers.csv", help="name of output csv table")

args = parser.parse_args()
species_collection = pd.read_csv(args.t, index_col="Abbreviation")
base_directory = args.d
proteins_ordered = args.p.split(",")
output = args.o

# Select the most recent orths file in each protein directory
texts = []
for root, dirs, files in os.walk(base_directory, topdown=False):
    if root == base_directory: 
        continue
    orth_files = [f for f in files if 'euk5_orths' in f and '.txt' in f]
    if len(orth_files) > 0:
        orth_file_latest = sorted(orth_files, reverse=True)[0]
        orth_file_latest_path = f"{root}/{orth_file_latest}"
        texts.append(orth_file_latest_path)
    else:
        print(f"{root} contains no orths text file")

# Create a dictionary of the proteins with their most recent orths txt files
protein_texts = {textfile.split(".")[-2]:textfile for textfile in texts} # File name convention: prefix - protname - suffix (.txt)

# For each protein in the ordered proteins list, check if they are represented by an orths file
# If so, collect the occurrences of this protein in the species of eukarya.v5 and append that to the species dataframe
for protein in proteins_ordered:
    if protein in protein_texts:
        textfile = protein_texts[protein]
        if os.path.exists(textfile):
            species_collection[protein] = ""
            infile = open(textfile, 'r')
            flines = infile.readlines()
            sequence_ids = [l.rstrip("\n") for l in flines]
            hits_species = {sp:[] for sp in species_collection.index.tolist()}
            for s in sequence_ids:
                try:
                    s_species = re.search(r'(\D{6})(\d{6})', s)[1]
                    hits_species[s_species].append(s)
                except:
                    print(f"Error: no species found for {s}")
            for k, v in hits_species.items():
                species_collection.loc[k, protein] = ";".join(v)
            infile.close()
        else:
            print(f"text file not found: {textfile}")
    else:
        print(f"protein has no textfile: {protein}")

# Print the dataframe to a csv 
species_collection.to_csv(output, index=True)
