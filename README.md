# AI-based Prediction of Protein Corona Composition on DNA Nanostructures
==================

Updated August 15, 2024

This repository contains the code and data necessary to reproduce figures from "AI-based Prediction of Protein Corona Composition on DNA Nanostructures" by Jared Huzar*, Roxana Coreas* Markita P. Landry, and Grigory Tikhomirov. This article is available in preprint form on BioxRiv at https://www.biorxiv.org/content/10.1101/2024.08.25.609594v1.

Briefly, we develop a machine learning model that can determine whether a protein will be present in the protein corona of different nanostructures based on the design features of the DNA nanostructures and the functional, structural, and sequence properties of the protein. We also identify the properties of the proteins and the nanostructures that govern their adsorption.

*These authors contributed equally.

## Running the code

To run the code, please clone this GitHub repository and from the main directory run the 'Simplified_Complete_Code.ipynb' file on Jupyter Notebook to do each analysis step-by-step including the construction of the database. Alternatively, the entire analysis can be run from command line by running  'Simplified_Complete_Code.py'. Results will be output in the 'Results' folder, and intermediate data files will be in the 'AnalyzedDataTest' folder.

The code was developed and tested on Windows 11 with python3.10. The following modules are required for running the entire analysis:
- Numpy
- Pandas
- Seaborn
- Math
- Quantiprot
- Matplotlib
- Stats
- XGBoost
- Biopython
- re
- requests
- sklearn
- shap
- scipy
- upsetplot
- matplotlib-venn


Individual functions can also be ran on their own on the command line from the main directory by running the intended python files found in either the 'AnalysisScripts' directory or the 'DatabaseFormation' directory.

## Scripts and data contained

The following subdirectories are housed within this repository:
- AnalysisScripts
	- Contains the scripts that produce the figures, run the machine learning model, and broadly analyze the nanostructures protein corona compositions.
- AnalyzedData
	- Contains the processed data which is required for the analysis
- Data
	- Contains the LCMS data, the Netsurf data, and the origami features database
- DatabaseFormation
	- Contains the python scripts that produce the database of protein properties and merges that data with the LCMS and origami data.
- Results
	- Contains figures and tables of statistics produced by the scripts in the AnalysisScripts directory.


In the subdirectories, we have placed a readme file explaining the data, results, and/or python scripts within the respective directories.

