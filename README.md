# SMCevolution
**Evolutionary histories and current-day distributions of eukaryotic SMC complexes' proteins**

## proteins
* Lists (\*.txt) and FASTA files (\*.fa) of SMC complexes' proteins
  * Filename convention: [dataset]\_[protein set type].[protein name].[txt or fa]; e.g. euk5_orths2.SMC1.fa, 'euk5' indicating eukarya.v5
  * Protein set type: 'orths' (sequences belonging to the orthogroup within eukaryotes) or 'homs' (broader set of homologs, e.g. also containing eukaryotic outparalogs or prokaryotic homologs)
* HMM profiles used to establish the orthogroup or to gather homologs for the phylogeny
* Phylogenies used to establish the orthogroup, including:
  * Treefile(s)
  * Logfile
  * Multiple sequence alignment

## protein_families
* SMC, kleisin and kite: Multiple sequence alignments (unfiltered and filtered, the latter were used for the phylogenies) and phylogenies of protein families studied to examine the origins of the eukaryotic SMC complexes, including constraint phylogenies used for topology testing
* Hawks: 
  * Curated profile HMM files used in profile-versus-profile homology searches with HHsearch, including also MAU2, WAPL, Nse5 and Nse6 (other alpha-solenoid domain proteins)
  * PyMOL session file with alignments of predicted hawk structures of *Homo sapiens* and *Arabidopsis thaliana*, derived from the AlphaFold repository (alphafold.ebi.ac.uk/, last accessed 23 December 2023)
* Nse56: 
  * AlphaFold2-predicted structures for a subset of Nse5 and Nse6 orthologs, used as input for structural alignments; both provided in full-length (separate Nse5 and Nse6 zipped archives contain full-length predictions, the combined folder contains the cut alpha-solenoid domains)
  * Structure alignment created with TMalign, respresenting the best Nse5-Nse6 alignment across all Nse5-Nse6 pairs of AlphaFold2 predicted structures

## euk5proteomes
* Tab-separated files (\.csv) listing the eukarya.v5 ('euk5') proteome dataset
* Species phylogenies (\.nwk and \.nexus) of the taxa eukarya.v5 - so far topologies only

## phylogenetic_profiles
* Tab-separated files (\.csv) showing the number of orthologs (members of the orthogroups) in each taxon of eukarya.v5

## py_scripts
* General-usage python3 scripts, e.g. to generate phylogenetic profile tables

## gtdb_selection
* Metadata for the included assemblies from GTDB (archaea, bacteria)

