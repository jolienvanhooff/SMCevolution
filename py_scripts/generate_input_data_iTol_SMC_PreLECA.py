#!/usr/bin/env python3

"""
generate_input_data_iTol_SMC_PreLECA.py

This script generates visualization datasets for the Interactive Tree of Life (iTOL) to annotate 
phylogenetic trees of SMC complex subunits (SMC, kleisin and kite proteins). It focuses on the evolutionary 
relationships of SMC complex subunits across the three domains of life (Bacteria, Archaea, Eukaryota), with 
particular emphasis on the eukaryotic SMC complexes (condensin I, condensin II, cohesin, SMC5/6).

The script processes a phylogenetic tree file and creates several iTOL annotation files:
1. Branch coloring based on taxonomic domain
2. Relabeling of taxa with relevant taxonomic and protein information
3. Highlighting of paralogs in the tree
4. Annotation of SMC complex memberships for eukaryotic proteins

This analysis is particularly useful for studying the evolution of SMC complex subunits before LECA 
(Last Eukaryotic Common Ancestor), helping to understand the origin and diversification of 
these essential chromosome organization proteins.

Required input:
- Phylogenetic tree file in Newick format
- Metadata tables for Archaea, Bacteria, and Eukaryota
- Optional text files with protein family memberships

Output:
- Reformatted tree file
- Multiple iTOL annotation datasets for visualization

Example usage: generate_input_data_iTol_SMC_PreLECA.py -t ../protein_families/Kite/euk5_homs6.Nse1_Nse3.ginsi_dash.cut.gappyout.drop.iqtree.treefile -r Arch_Asgardarchaeota_Heimdallarchaeia_GB_GCA_001940645.1_MDVS01000047.1_5 Arch_Altiarchaeota_Altiarchaeia_GB_GCA_016935655.1_JAFGQM010000011.1_41 -p clade.Nse1.txt clade.Nse3.txt clade.Kite_unknown_archaea.txt clade.Nse1_Nse3_related.txt  clade.Nse3_related.txt
The clade files contains identifiers of proteins belonging to a particular subfamilies, or clade, in the phylogeny, either prokaryotic or eukaryotic. For instance, clade.Nse1.txt contains a list with the identifiers of the eukaryotic Nse1 proteins: 
HALSEO004537_Nse1
ROZALL004763_Nse1
TRIMAR002709_Nse1
(list continues with all Nse1 orthologs in this phylogeny)
"""


from ete3 import Tree
import argparse
import pandas as pd
import os
from functions import checktrailingslash, make_list_from_lines, find_largest_common_prefix
import re
import collections


# define the color scheme for each leaf prefix
color_scheme = {
    "Bacteria": "#0044AD",
    "Archaea": "#7F6000",
    "Eukaryota": "#C90067"
}

# define protein complex memberships for SMC complexes (eukaryotes)
complex_members = {
    "SMC2": "Condensin", 
    "SMC4": "Condensin",
    "SMC3": "Cohesin",
    "SMC1": "Cohesin",
    "SMC5": "SMC5/6",
    "SMC6": "SMC5/6",
    "CAPH": "CondensinI",
    "CAPH2": "CondensinII",
    "Scc1": "Cohesin",
    "Scc1_Rec8": "Cohesin",
    "Rec8": "Cohesin",
    "Nse4": "SMC5/6",
    "Nse1": "SMC5/6",
    "Nse3": "SMC5/6"
}

# define a separate colour scheme for SMC complex protein membership (eukaryotes)
colors_complex_members = {
    "SMC2": "#ddcc77", 
    "SMC4": "#ddcc77",
    "SMC3": "#4477aa",
    "SMC1": "#4477aa",
    "SMC5": "#aa3377",
    "SMC6": "#aa3377",
    "CAPH": "#999933",
    "CAPH2": "#117733",
    "Scc1": "#4477aa",
    "Scc1_Rec8": "#4477aa",
    "Rec8": "#4477aa",
    "Nse4": "#aa3377",
    "Nse1": "#aa3377",
    "Nse3": "#aa3377"
}

def simplify_leaf_names(tree):
    for leaf in tree:
        leaf.name = leaf.name.split("/")[0]


def reroot_tree(tree, root_leaves):
    if len(root_leaves) == 1:
        outgroup_node = tree.search_nodes(name=root_leaves[0])[0]
        tree.set_outgroup(outgroup_node)
    elif len(root_leaves) == 2:
        ancestor = tree.get_common_ancestor(root_leaves)
        tree.set_outgroup(ancestor)


def clean_tree(tree):
    """Remove branches corresponding to an AlphaFold predicted structure from the phylogeny"""
    if any(l.name.startswith("AF-") for l in tree.get_leaves()):
        leaves_to_preserve = [l for l in tree.get_leaves() if l.name.startswith("AF-") == False]
        tree.prune(leaves_to_preserve, preserve_branch_length=True)


def get_furthest_descendants(node):
    if node.is_leaf() == False:
        child1, child2 = node.get_children()[0], node.get_children()[1]
        # get a random leaf of each child
        leaf1 = child1.get_leaf_names()[0]
        leaf2 = child2.get_leaf_names()[-1]
        return [leaf1, leaf2]
    else:
        return [node.name]


def get_protein_memberships_from_txtfiles(filelist):
    seqname_to_protgroup = {}
    try: 
        for f in filelist:
            prot_name = f.split(".")[-2] # filenames end with '.txt'
            members = make_list_from_lines(f)
            for member in members:
                seqname_to_protgroup[member] = prot_name
    except TypeError: 
        print("No txt files with protein memberships provided")
    return seqname_to_protgroup


def label_leaves(tree, arch_metadata, bact_metadata, euk_metadata, protein_memberships):
    """Assign labels to leaves (protein name (of group), domain, phylum, class and taxonomy for prokaryotes and the protein for eukaryotes)"""
    for l in tree.get_leaves():
        l.protein = protein_memberships[l.name] if l.name in protein_memberships else ""

        if l.name.startswith("Arch"):
            # Annotate archaea
            l.domain = "Archaea"
            l.phylum, l.clas = l.name.split("_")[1], l.name.split("_")[2]
            l.accession = re.search(r'_((GB|RS)_[^.]+\.\d+)_', l.name).group(1)
            if l.accession not in arch_metadata.index:
                print(f"Error: {l.accession} not found in metadata Archaea")
            else: 
                l.taxonomy = arch_metadata.loc[l.accession, "gtdb_taxonomy"]
                l.species = l.taxonomy.split("__")[-1]

        elif l.name.startswith("Bact"):
            # Annotate bacteria
            l.domain = "Bacteria"
            l.phylum, l.clas = l.name.split("_")[1], l.name.split("_")[2]
            l.accession = re.search(r'_((GB|RS)_[^.]+\.\d+)_', l.name).group(1)
            if l.accession not in bact_metadata.index:
                print(f"Error: {l.accession} not found in metadata Bacteria")
            else:
                l.taxonomy = bact_metadata.loc[l.accession, "gtdb_taxonomy"]
                l.species = l.taxonomy.split("__")[-1]
        
        else:
            # Annotate eukaryotes; replace eukaryotic label (is behind the protein identifier)
            l.domain = "Eukaryota"
            if re.search(r'^(\D{6})(\d{6})', l.name):
                l.protein = l.name.split("_")[1]
                l_species_acronym = re.search(r'^(\D{6})(\d{6})', l.name).group(1)
                l.species = euk_metadata.loc[l_species_acronym, "Scientific name"]
                l.rel_clade = euk_metadata.loc[l_species_acronym, "relevant taxonomy"]


def label_internal_nodes(tree):
    """Assign labels to internal nodes: protein names corresponding to the (orthologous) groups and the clade name and rank for prokaryotic clades only"""
    hierarchy = collections.OrderedDict({'s':'species', 'g':'genus', 'f':'family', 'o':'order', 'c':'class', 'p':'phylum', 'd':'domain'})
    for n in tree.iter_descendants("preorder"):
        if n.is_leaf() == False:
            leaves = n.get_leaves()
            # First get protein name; either use the protein that all leaves are labelled with, or the one dat some leaves are labelled with, whereas others are empty
            leaves_proteins = [l.protein for l in leaves]
            if all(prot == leaves_proteins[0] for prot in leaves_proteins):
                n.protein = leaves_proteins[0]
            # Make an exception for Scc1 and Rec8: known to be interspersed
            elif set(leaves_proteins) == {"Scc1", "Rec8"}:
                n.protein = "Scc1_Rec8"
            elif len(set(leaves_proteins))==2 and "" in set(leaves_proteins) and any(s != "" for s in set(leaves_proteins)):
                n.protein = next(s for s in set(leaves_proteins) if s != "")
            else: 
                n.protein = ""
            # Then get the domain, an the clade and rank according to the leaf taxonomies (prokaryotic only)
            leaves_domain = [l.domain for l in leaves]
            if all(item == leaves_domain[0] for item in leaves_domain):
                n.domain = leaves_domain[0]
                if n.domain == "Archaea" or n.domain == "Bacteria":
                    leaves_taxonomies = [l.taxonomy for l in leaves]
                    n_taxonomy_full = find_largest_common_prefix(leaves_taxonomies)
                    ## If not shared up to the species level, the pattern will end with the rank and "__", which should be removed
                    if n_taxonomy_full.endswith("__"):
                        n_taxonomy_full = n_taxonomy_full[0:-4]
                    # Select the lowest level
                    n_taxonomy_specific = n_taxonomy_full.split(";")[-1]
                    n.rank, n.clade = n_taxonomy_specific.split("__")
                    n.rank = hierarchy[n.rank]
                else: 
                    n.clade = ""
                    n.rank = ""
            else:
                n.domain = ""
                n.clade = ""
                n.rank = ""
    # Iterate once more to label (protein name, species domain) of unannotated nodes based on parent and tips (how common? probably not very given the strict requirements of the internal node labelling)
    for n in tree.iter_descendants("preorder"):
        if n.is_leaf() == False:
            if n.protein == "":
                try:
                    p_protein = n.up.protein
                    l_proteins = [c.protein for c in n.get_leaves()]
                    # Check if there's one match between parent and children
                    if any(l_protein == p_protein for l_protein in l_proteins):
                        n.protein = p_protein
                except AttributeError:
                    ""
            if n.domain == "":
                try: 
                    p_domain = n.up.domain
                    l_domains = [c.domain for c in n.get_leaves()]
                    # Check if there's one match between parent and children
                    if any(l_domain == p_domain for l_domain in l_domains):
                        n.domain = p_domain
                except AttributeError:
                    pass
                

def generate_itol_dataset_branch_colours(tree):
    """Generate iTol dataset from a tree and a root prefix"""
    # create the iTol dataset in comma-separated format
    dataset = []
    # Find monophyletic groups of bacteria, archaea and eukaryotes and colour them (the ancestor and the entire clade emanating from it)
    for dom in color_scheme.keys():
        for n in tree.get_monophyletic(values=[dom], target_attr="domain"):
            descendants = get_furthest_descendants(n)
            dataset.append(f"{'|'.join(descendants)},branch,clade,{color_scheme[dom]},1,normal")
    # Also colour those internal nodes that aren't monophyletic with regard to the domain, but that do have an annotated domain (test)
    for n in tree.iter_descendants("preorder"):
        if n.is_leaf() == False:
            for dom in color_scheme.keys():
                if n.domain == dom:
                    if n.check_monophyly(values=[dom], target_attr="domain") == False:
                        descendants = get_furthest_descendants(n)
                        dataset.append(f"{'|'.join(descendants)},branch,clade,{color_scheme[dom]},1,normal")
    return dataset


def generate_itol_dataset_new_names(tree):
    """Generate iTol dataset that renames the taxa"""
    dataset = []
    for l in tree.get_leaves():
        if l.name.startswith(("Bact", "Arch")):
            dataset.append(f"{l.name},{l.phylum}_{l.species.replace(' ', '_')}")
        else: 
            dataset.append(f"{l.name},{l.rel_clade}_{l.species.replace(' ', '_')}")
    # Label the internal nodes with the protein name (all domains) and/or the clade (prokaryotes only) - depending on what's available
    for n in tree.iter_descendants("preorder"):
        if n.is_leaf() == False:
            descendants = get_furthest_descendants(n)
            if n.protein != "" and n.clade != "":
                dataset.append(f"{'|'.join(descendants)},{n.clade}_{n.protein}")
            elif n.protein != "":
                dataset.append(f"{'|'.join(descendants)},{n.protein}")
            elif n.clade != "":
                dataset.append(f"{'|'.join(descendants)},{n.clade}")
    return dataset


def get_color_list(num_colors):
    import colorsys
    # Generate equally spaced hues
    hues = [i / num_colors for i in range(num_colors)]
    # Convert hues to RGB values
    rgb_values = [colorsys.hsv_to_rgb(hue, 1, 1) for hue in hues]
    # Convert RGB values to hex codes
    hex_list = ['#%02x%02x%02x' % tuple(int(round(255*val)) for val in color) for color in rgb_values]
    return hex_list


def generate_itol_dataset_paralogs(tree, arch_metadata, bact_metadata):
    """Add coloured symbols to paralog leaves"""
    dataset = []
    paralog_taxa = []
    # First get the taxa that have paralogs in the tree
    for l in tree.get_leaves():
        if l.domain == "Archaea" and l.accession in arch_metadata.index:
            if len(tree.search_nodes(domain="Archaea", accession = l.accession)) > 1 and l.name not in paralog_taxa:
                paralog_taxa.append(l.accession)
        elif l.domain == "Bacteria" and l.accession in bact_metadata.index:
            if len(tree.search_nodes(domain="Bacteria", accession = l.accession)) > 1 and l.name not in paralog_taxa:
                paralog_taxa.append(l.accession)
    paralog_taxa = list(set(paralog_taxa))
    paralog_taxa.sort()
    # Add star symbols for each taxon with duplicates
    hex_range = get_color_list(len(paralog_taxa))
    for i, taxon in enumerate(paralog_taxa):
        taxon_nodes = tree.search_nodes(accession = taxon)
        for paralog in taxon_nodes:
            dataset.append(f"{paralog.name},3,1,{hex_range[i]},1,0.8")
    return(dataset)


def generate_itol_dataset_membership(tree):
    """Collect SMC complex memberships for leaves and for internal nodes; create vertically coloured bars - add complex for internal nodes"""
    dataset = []
    for l in tree.get_leaves():
        if l.protein != "" and l.protein in colors_complex_members:
            prot_complex_color = colors_complex_members[l.protein]
            dataset.append(f"{l.name},{prot_complex_color}")
    for n in tree.iter_descendants("preorder"):
        if n.is_leaf() == False:
            if n.protein != "" and n.protein in colors_complex_members:
                descendants = get_furthest_descendants(n)
                prot_complex_color = colors_complex_members[n.protein]
                prot_complex = complex_members[n.protein]
                dataset.append(f"{'|'.join(descendants)},{prot_complex_color},{prot_complex}")
    return dataset


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate iTol datasets to annotate a gene phylogeny')
    parser.add_argument('-t', metavar='tree_path', type=str, help='Path to the input newick tree file')
    parser.add_argument('-ma', metavar='archaea', type=str, help='Path to metadata table of archaeal lineages', default="../gtdb_selection/ar53_metadata_r207.qscore.family_representative.csv")
    parser.add_argument('-mb', metavar='bacteria', type=str, help='Path to metadata table of bacterial lineages', default="../gtdb_selection/bac120_metadata_r207.qscore.family_representative.csv")
    parser.add_argument('-me', metavar='eukaryota', type=str, help='Path to metadata table for eukaryotes', default="../euk5proteomes/Euk5FinalSet.adjust.busco.euk5_tree_abbrev.csv")
    parser.add_argument('-o', metavar='output_dir', type=str, help='output directory - current working directory if not specified')
    parser.add_argument('-r', metavar='root_leaves', nargs='+', type=str, help='List of leaf names for rooting the tree')
    parser.add_argument('-p', metavar='protein_membership', nargs='+', type=str, help='List of text files containing subfamily memberships of (prokaryotic) sequences - suffix should be ".txt"')
    args = parser.parse_args()

    # Get the output directory
    if args.o == None:
        outdir = checktrailingslash(os.getcwd())
    else:
        outdir = checktrailingslash(args.o)

    # Get the basename of the input tree
    tree_basename = os.path.basename(args.t)
    
    # Load the tree from file
    tree = Tree(args.t)

    # Simplify leaf names
    simplify_leaf_names(tree)

    # Reroot the tree
    reroot_tree(tree, args.r)

    # Remove sequences from AlphaFold structures from the tree
    clean_tree(tree)

    # Write new tree to newick
    tree.write(outfile=f"{outdir}{tree_basename}.reformatted", format=0)
    
    # Load species metadata
    archaea_metadata = pd.read_csv(args.ma, index_col="accession")
    bacteria_metadata = pd.read_csv(args.mb, index_col="accession")
    eukaryota_metadata = pd.read_csv(args.me, index_col="Abbreviation")
    
    # Load protein memberships (for prokaryotes)
    protein_memberships_from_txt = get_protein_memberships_from_txtfiles(args.p)
    
    # Label the leaves according to their domain and lower taxonomy or protein
    label_leaves(tree, archaea_metadata, bacteria_metadata, eukaryota_metadata, protein_memberships_from_txt)

    # Label the internal nodes: protein name and taxonomy (the latter for prokaryotes only)
    label_internal_nodes(tree)

    # Generate and write the iTol dataset to branch colours according to the species domain (Eukaryota, Archaea, Bacteria)
    dataset_branch_colours = generate_itol_dataset_branch_colours(tree)
    with open(f"{outdir}{tree_basename}.reformatted.iTOL_domain.dataset.txt", "w") as file:
        file.write("DATASET_STYLE\n")
        file.write("SEPARATOR COMMA\n")
        file.write("DATASET_LABEL,Domain\n")
        file.write("COLOR,#ffff00\n")
        file.write("LEGEND_TITLE,Domain\n")
        file.write("LEGEND_POSITION_X,100\n")
        file.write("LEGEND_POSITION_Y,100\n")
        file.write("LEGEND_HORIZONTAL,0\n")
        file.write("LEGEND_SHAPES,1,1,1\n")
        file.write("LEGEND_COLORS,#C90067,#7F6000,#0044AD\n")
        file.write("LEGEND_LABELS,Eukaryota,Archaea,Bacteria\n")
        file.write("LEGEND_SHAPE_SCALES,1,1,1\n")
        file.write("DATA\n")
        for item in dataset_branch_colours:
            file.write(f"{item}\n")    

    # Generate and write the iTol dataset for new names
    dataset_labels = generate_itol_dataset_new_names(tree)
    with open(f"{outdir}{tree_basename}.reformatted.iTOL_labels.dataset.txt", "w") as file:
        file.write("LABELS\n")
        file.write("SEPARATOR COMMA\n")
        file.write("DATA\n")
        for item in dataset_labels:
            file.write(f"{item}\n")
    
    # Generate and write the iTol dataset for paralogs 
    dataset_paralogs = generate_itol_dataset_paralogs(tree, archaea_metadata, bacteria_metadata)
    with open(f"{outdir}{tree_basename}.reformatted.iTOL_paralogshapes.dataset.txt", "w") as file:
        file.write("DATASET_SYMBOL\n")
        file.write("SEPARATOR COMMA\n")
        file.write("DATASET_LABEL,Paralogs\n")
        file.write("COLOR,#AC3A6D\n")
        file.write("LEGEND_TITLE,Paralogs\nLEGEND_POSITION_X,80\nLEGEND_POSITION_Y,80\nLEGEND_HORIZONTAL,0\nLEGEND_SHAPES,3\nLEGEND_COLORS,#AC3A6D\nLEGEND_LABELS,paralog\nLEGEND_SHAPE_SCALES,1\nLEGEND_SHAPE_INVERT,0\n")
        file.write("MAXIMUM_SIZE,10\n")
        file.write("#GRADIENT_FILL,1\n")
        file.write("DATA\n")
        for item in dataset_paralogs:
            file.write(f"{item}\n")
    
    # Generate and write the iTOL dataset for SMC complex memberships - eukaryotic proteins only - create a colored strip dataset for visualizing these memberships
    dataset_complex_membership = generate_itol_dataset_membership(tree)
    with open(f"{outdir}{tree_basename}.reformatted.smc_complex_memberships.txt", "w") as file:
        file.write(f"DATASET_COLORSTRIP\n")
        file.write(f"SEPARATOR COMMA\n")
        file.write(f"DATASET_LABEL,SMC_complex\n")
        file.write(f"COLOR,#8fce00\n")
        file.write(f"COLOR_BRANCHES,0\n")
        file.write(f"LEGEND_TITLE,SMC_complexes\n")
        file.write(f"LEGEND_POSITION_X,60\n")
        file.write(f"LEGEND_POSITION_Y,60\n")
        file.write(f"LEGEND_HORIZONTAL,0\n") 
        file.write(f"LEGEND_SHAPES,1,1,1,1,1\n")
        file.write(f"LEGEND_COLORS,{colors_complex_members['SMC2']},{colors_complex_members['CAPH']},{colors_complex_members['CAPH2']},{colors_complex_members['SMC1']},{colors_complex_members['SMC5']}\n") 
        file.write(f"LEGEND_LABELS,Condensin,CondensinI,CondensinII,Cohesin,SMC5/6\n") 
        file.write(f"LEGEND_SHAPE_SCALES,1,1,1,1,1\n") 
        file.write(f"STRIP_WIDTH,25\n")
        file.write(f"MARGIN,0\n")
        file.write(f"BORDER_WIDTH,1\n")
        file.write(f"BORDER_COLOR,#000000\n")
        file.write(f"COMPLETE_BORDER,1\n")
        file.write(f"SHOW_INTERNAL,0\n")
        file.write(f"SHOW_STRIP_LABELS,1\n")
        file.write(f"STRIP_LABEL_POSITION,center\n")
        file.write(f"STRIP_LABEL_SIZE_FACTOR,0.5\n")
        file.write(f"STRIP_LABEL_ROTATION,0\n")
        file.write(f"STRIP_LABEL_SHIFT,0\n")
        file.write(f"STRIP_LABEL_COLOR,#000000\n")
        file.write(f"SHOW_LABELS,0\n")
        file.write(f"DATA\n")
        for item in dataset_complex_membership:
            file.write(f"{item}\n")

                   













