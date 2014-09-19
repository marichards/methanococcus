function model = addKbaseAliases(model)

%Take in the M.maripaludis model and add the proper aliases

%Load the compiled dictionaries
load('compound_reaction_dicts.mat')

%Now compounds and reactions are in, both structs



%COMPOUNDS
%compounds has 3 fields: seedIDs, names, keggIDs

%Add a field for SEED IDs...model.metSEEDID
%Compound IDs for Kegg go in model.metKEGGID


%REACTIONS
%reactions has 4 fields: seedIDs, names, ECs, keggIDs

%Add a field for reaction kegg IDs...model.rxnKEGGID
%EC numbers go in model.rxnECNumbers