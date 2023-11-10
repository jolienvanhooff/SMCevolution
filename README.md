# SMCevolution
**Evolutionary histories and current-day distributions of eukaryotic SMC complexes' proteins**

***Warning: current version contains outdated orthologous groups for Sororin, Nse5 and Nse6 (groups for eukarya.v5.1 instead of the most recent eukarya.v5.2)***

## proteins
* Lists (\*.txt) and FASTA files (\*.fa) of SMC complexes' proteins
  * Filename convention: [dataset]\_[protein set type].[protein name].[txt or fa]; e.g. euk5_orths2.SMC1.fa, 'euk5' indicating eukarya.v5
  * protein set type: 'orths' (sequences belonging to the orthogroup within eukaryotes) or 'homs' (broader set of homologs, e.g. also containing eukaryotic outparalogs or prokaryotic homologs)
* HMM profiles used to establish the orthogroup or to gather homologs for the phylogeny
* Phylogenies used to establish the orthogroup, including:
  * (annotated) treefile
  * logfile
  * multiple sequence alignment

## protein_families
* MSAs and phylogenies of protein families studied to examine the origins of the eukaryotic SMC complexes

## alveolata
* Subproject: contains data similar to the main directory, most importantly 'proteins' based on searches across Jolien's SAR dataset (418 taxa)

## euk5proteomes
* Tab-separated files (\.csv) listing the eukarya.v5 ('euk5') proteome dataset
* Species phylogenies (\.nwk and \.nexus) of the taxa eukarya.v5 - so far topologies only

## phylogenetic_profiles
* Tab-separated files (\.csv) showing the number of orthologs (members of the orthogroups) in each taxon of eukarya.v5
* Species phylogeny figures showing the distributions of SMC complexes' proteins across eukarya.v5

## py_scripts
* General-usage python3 scripts, e.g. to generate phylogenetic profile tables

## gtdb_selection
* Metadata for the included assemblies from GTDB (archaea, bacteria)

