Contents of the Model Folder:

OLD MODEL VERSIONS(Directory)
Contains the original SBML file and some early model manipulations

original_model.mat
Contains the original model received from Matt Benedict

{date}_model.mat
Contains the current version of the model

{date}_dde_model.mat
Contains the current version of the model with double dead ends removed using removeDoubleDeadEnds(model) code.  Also contains the removed double dead end metabolites and reactions

{date}_sde_model.mat 
Contains the current version of th emodel with double AND single dead ends removed using removeDeadEnds(dde_model) code.  Also contains removed single dead end metabolites and reactions

All Model Events.txt
Human-readable file that chronicles all the changes, in order, that we've made to the original model.  Should align with createLatestModel.m

createLatestModel.m
Script that creates the current model from the original.  Should be changed and commented when we want to make changes to the model