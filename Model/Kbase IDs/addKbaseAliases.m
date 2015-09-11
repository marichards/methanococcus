function model = addKbaseAliases(model)

%Take in the M.maripaludis model and add the proper aliases

%Load the compiled dictionaries
load('compound_reaction_dicts.mat')

%Now compounds and reactions are in, both structs

%COMPOUNDS
%compounds has 5 fields: seedIDs, names, keggIDs,names_c0,names_e0

%Add a field for SEED IDs...model.metSEEDID
[m,n]=size(model.mets);
model.metSEEDID = cell(m,n);
%Compound IDs for Kegg go in model.metKEGGID
%Do it for c0 ones first:
[~,model_idx,dict_idx] = intersect(model.mets,compounds.names_c0);
for i=1:length(model_idx)
    %For each compound, assign to the MODEL index the values from the DICT
    %index
    model.metSEEDID{model_idx(i)}=compounds.seedIDs{dict_idx(i)};
    model.metKEGGID{model_idx(i)}=compounds.keggIDs{dict_idx(i)};
end

%Do the same for e0 ones
[~,model_idx,dict_idx] = intersect(model.mets,compounds.names_e0);
for i=1:length(model_idx)
    %For each compound, assign to the MODEL index the values from the DICT
    %index
    model.metSEEDID{model_idx(i)}=compounds.seedIDs{dict_idx(i)};
    model.metKEGGID{model_idx(i)}=compounds.keggIDs{dict_idx(i)};
end

%REACTIONS
%reactions has 5 fields: seedIDs, names, ECs, keggIDs, seedIDs_w_c0
%The last one is our comparison and should have 686 in common
%I added c0 to it!

%%%IT WORKS TO COMPARE! HUZZAH!
[~,model_idx,dict_idx] = intersect(model.rxns,reactions.seedIDs_w_c0);

%Now I have the list of reactions that are in the model AND the dictionary,
%plus the index for each one of them, though I'm only keeping the indices

%Add the info:
%Add a field for reaction kegg IDs...model.rxnKEGGID
[m,n]=size(model.rxns);
model.rxnKEGGID = cell(m,n);
%EC numbers go in model.rxnECNumbers

%Loop through:
for i=1:length(model_idx)
    %For each reaction, assign to the MODEL index the values from the DICT
    %index
    model.rxnKEGGID{model_idx(i)}=reactions.keggIDs{dict_idx(i)};
    model.rxnECNumbers{model_idx(i)}=reactions.ECs{dict_idx(i)};    
end

% We missed some Seed IDs for compounds and subsystems for reactions
load('metIDs.mat') 
[~,idxA,idxB] = intersect(model.mets,metIDs.name);
model.metSEEDID(idxA)=metIDs.ID(idxB);

% Add subsystems in similar fashion
load('subsystems.mat')
[~,idxA,idxB] = intersect(model.rxns,subsystems.IDs);
model.subSystems(idxA)=subsystems.subs(idxB);

