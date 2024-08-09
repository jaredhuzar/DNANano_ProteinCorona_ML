# AI-based Prediction of Protein Corona Composition on DNA Nanostructures
==================
Updated August 9, 2024

## AnalysisScripts

This directory contains scripts that create, test, and evaluate the machine learning model. In addition, they perform several other analyses elucidating different features of DNA nanostructures' coronas.

- RunModelAnalysisPresent.py
	- Develops an XGBoost machine learning model that predicts whether or not a protein will be adsorbed in the corona of different DNA nanostructures. Tests and evaluates this model.

- RunModelEnriched.py
	- Develops an XGBoost machine learning model that predicts whether or not a protein will be enriched in the corona of different DNA nanostructures. Tests and evaluates this model.

- PerformanceAcrossSubsets.py
	- Independently trains and evaluates XGBoost model on only DNA nanostructure data, only protein functional data, and only protein structural data and evaluates the effectiveness of each of these data subsets in predicting whether a protein will be present in nanostructures' coronas.

- Modified_vsUnmodified.py
	-Compares the protein composition between DNA only nanostructures, and those modified with either a polymer coating or cholesterol.

- FunctionalEnrichment.py
	- Evaluates whether proteins with different molecular functions are more present in the coronas than expected by random chance.

- EndocyticProteinEnrichment.py
	- Evaluates the enrichment of proteins associated with endocytosis in the corona of different nanostructures.

- DifferentProteinProperties.py
	- Evaluates the difference between proteins universally present in all nanostructure coronas and those not found in any nanostructure coronas.

- CommonProteins.py
	- Evaluates the compositional similarity between coronas of different nanostructures.

- CoatedvsUncoatedModels.y
	- Trains and evaluates models for predicting protein in-corona presence on coated and bare structures separately.

- ModelComparison.py
	- Evaluates the performance of three different ML models for the classification of in-corona presence.

- RunClassifier_EachNanostructure
	- Trains and evaluates XGBoost models trained on each nanostructure individually.
