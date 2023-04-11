# SMCevolution
**Evolutionary histories and current-day distributions of eukaryotic SMC complexes' proteins**

***Warning: current version contains outdated orthologous groups for Sororin, Nse5 and Nse6 (groups for eukarya.v5.1 instead of the most recent eukarya.v5.2)***

*NB: typically, this repository serves to share and visualize results and, importantly, to keep track of the versions of such results. It is not intended to use as a location for storing all sorts of intermediate files*

## proteins
* Lists (\*.txt) and FASTA files (\*.fa) of SMC complexes' proteins
  * Filename convention: [dataset]\_[protein set type].[protein name].[txt or fa]; e.g. euk5_orths2.SMC1.fa, 'euk5' indicating eukarya.v5
  * protein set type: 'orths' (sequences belonging to the orthogroup within eukaryotes) or 'homs' (broader set of homologs, e.g. also containing eukaryotic outparalogs or prokaryotic homologs)
* HMM profiles used to establish the orthogroup or to gather homologs for the phylogeny
* Phylogenies used to establish the orthogroup, including:
  * (annotated) treefile
  * logfile
  * multiple sequence alignment
* *other data, e.g. for pre-LECA examinations - maybe it is most practical if we add a subdirectory for these*

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

## co-evolution *to be added*

## .gitignore
* Add anything to be ignored in version controlling

## External links / other platforms
* [Google Drive directory](https://drive.google.com/drive/folders/10_zQJfThDdbN8nEHA8leA_tiyxAh5Fh-)
  * [Manuscript_outline](https://docs.google.com/document/d/1BkOMaUu7r-3rs05RzT2_LqnterpXH63XLRaZpAGHmeM/edit)
  * [General notes/observations](https://docs.google.com/document/d/1uRSr-7Q_5-_9Sp_bZsohDjWt9PrOxKg1IUqID9gxZLE/edit)
* [Collected literature](https://paperpile.com/shared/Uow2va)
* [Slack channel](https://deemteamworkspace.slack.com/archives/C04579BJXR8)
* [Google Spaces](https://mail.google.com/chat/u/0/#chat/space/AAAAiR1Ycrk) *obsolete*
