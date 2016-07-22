## Model

This directory contains the most recent iteration of the *M. maripaludis* model, marked with the date it was created. It also contains the exact iMR539 model from publication in SBML and Matlab formats, the original automated model created using tools in Kbase, and a directory containing all versions of the model that came in between. 

### iMR539.xml

Contains the iMR539 model in SBML, compatible with all standard constraint based modeling platforms. 

### iMR539.mat

Contains the iMR539 model in a Matlab data structure. Notably, this model version contains additional metadata, including more alternate reaction and metabolite identifiers. It also includes the "freeEnergy" array, which catalogues standard free energies for exchange metabolites and is required for the "optimizeThermoModel.m" script and all functions that depend upon it. 

### original_model.mat

Contains the original draft reconstruction, created using automated reconstruction and likelihood-based gap-filling in Kbase. It lacks any manual curation, save for the HdrABC and Eha reactions added to the draft model. This model is essential for proper function of the createLatestModel.m script. 

### createLatestModel.m

This script takes in the original model from Kbase (original_model.mat) and modifies it to create the most current version of the model. Due to the living nature of the model, this script is constantly subject to augmentation and serves as a written record of every change made to the model. 

### {date}_model.mat

Contains the most current version of the model, tagged with the date that model was created. For the manuscript of iMR539, the equivalent model version is "2016_06_08_model.mat" (see Archived Models)

### Archived Models 

This directory contains all models that fall between the original_model.mat version and the most current version. All archived models are tagged with their date of creation. Notably, due to the changes in models over time, some of these models may not function correctly with the scripts in the "Matlab Scripts" directory. However, each model is still usable in conjunction with the COBRA Toolbox. 