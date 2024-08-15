# AI-based Prediction of Protein Corona Composition on DNA Nanostructures
==================

Updated August 15, 2024

## Database Formation

This directory contains the following scripts aimed at generating a database encompassing protein properties and nanostructure properties likely to drive protein adsorption to the nanoparticles. In addition, it merges this data with the experimentally measured LCMS protein data to create a dataset suitable for machine learning and other analyses.


- LCMS_BinaryExtraction.py
	- Creates dataframes that identify which proteins are present and enriched in the corona of each nanostructure.

- NetSurfAnalysis.py
	- Analyzes the data produced by NetSurf2.0

- UniProtID_Identification.py
	- Identifies the UniProt IDs associated with the proteins found in the serum and in the coronas.

- UniprotQuantiprot_Structural_Functional_Info.py
	- Scrapes UniProt and uses QuantiProt to quantify different structural and functional properties of the proteins.

- Merge_Protein_Data.py
	- Merges the protein properties derived from the different, aforementioned sources.

- Merge_Nano_Protein_LCMS.py
	- Merges protein properties, DNA nanostructure properties, and LCMS data into one dataset.

 

