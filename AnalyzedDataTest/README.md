# AI-based Prediction of Protein Corona Composition on DNA Nanostructures
==================

Updated August 15, 2024


## AnalyzedDataTest Folder

This repository contains the processed data files used to produce the analyses in the paper. This directory contains the following files:


- overallResEnriched8_6.csv
	- This csv file is necessary for running the model to predict whether proteins are enriched or depleted/absent on the protein corona of DNA nanostructures. The file contains DNA origami features, protein properties, and binary classifications of enriched or depleted/absent.

- overallResPresent8_6.csv
	- This csv file is necessary for running the model to predict whether proteins are present or absent on the protein corona of DNA nanostructures. The file contains DNA origami features, protein properties, and binary classifications of enriched or depleted/absent.

- ProStatistics8_7.csv
	- This csv file contains all of the properties for each protein.

- allProteins8_6.txt
	- This file contains all of the proteins analyzed and their respective amino acid sequences.

- flexibility.csv
	- This file contains the flexibility statistics for each protein.

- functionalInfo.csv
	- This file contains the information on each protein's functionality.

- FunctionList8_6.csv
	- This contains all of the functions considered in the analyses.

- InteractivityData.csv
	- A file containing the interactivity metric score for each protein

- ProData8_6.csv
	- File contains the protein information obtained from QuantiProt (CHECK)

- TranslatedPresent8_8.csv
	- Contains names of proteins used in the analyses and their respective UniProIDs

- aminoAcid.csv
	- Contains statistics on the amino acid composition of all analyzed proteins

- combos.csv
	- Contains statistics on the amino acid composition and their respective secondary structures for all analyzed proteins. (CHECK)

- exposedAminoAcid.csv
	- Contains statistics on the solvent-exposed amino acids for each protein.

- moreNetSurf.csv
	- Contains statistics derived from NetSurf 2.0 regarding structural properties of the proteins.

- secondaryAminoAcid.csv
	- Contains statistics on the amino acid composition and their respective secondary structures for all analyzed proteins.

- AllNetSurf.csv
	- File containing the combined data from 'moreNetSurf.csv', 'secondaryAminoAcid.csv', 'combos.csv', 'aminoAcid.csv', and 'exposedAminoAcid.csv'.

- EnrichedProteins8_4.csv
	- CSV file containing a list of all the proteins universally enriched.

- PresentProteins8_4.csv
	- A CSV file containing a list of all the proteins universally present.