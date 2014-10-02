function model = addModelInfo(model)

%Take in the other MM model, add various fields of info

%First run the code to add genes
model = addGenesToModel(model);

%Add bounds
model = addBounds(model);

%Remove that faulty gene
model = removeGene(model,'''''');

%Make biomass the objective
model = changeObjective(model,'R571',1);

%Read the rxn Info in
[~,rxns] = xlsread('supplementary 1.xls','Reactions','A2:A572');
[~,ECs] = xlsread('supplementary 1.xls','Reactions','H2:H572');
[~,rxnNames] = xlsread('supplementary 1.xls','Reactions','B2:B572');

%Read the metabolite Info in
[~,mets] = xlsread('supplementary 1.xls','Metabolites','A2:A606');
[~,metNames] = xlsread('supplementary 1.xls','Metabolites','B2:B606');
[~,metKEGGID] = xlsread('supplementary 1.xls','Metabolites','E2:E606');

%Intersect the metabolites, make the replacements

[~,model_idx,dict_idx] = intersect(model.mets,mets);
for i=1:length(model_idx)
    %For each compound, assign to the MODEL index the values from the DICT
    %index
    %This is for metNames and metKEGGID
    model.metNames{model_idx(i)}=metNames{dict_idx(i)};
    model.metKEGGID{model_idx(i)}=metKEGGID{dict_idx(i)};
end

%Intersect the reactions, make replacements
[~,model_idx,dict_idx] = intersect(model.rxns,rxns);

%Add a field for reaction EC numbers
[m,n]=size(model.rxns);
model.rxnECNumbers = cell(m,n);

for i=1:length(model_idx)
    %For each reaction, assign to the MODEL index the values from the DICT
    %index
    %This is for rxnNames and rxnECNumbers
    model.rxnNames{model_idx(i)}=rxnNames{dict_idx(i)};
    model.rxnECNumbers{model_idx(i)}=ECs{dict_idx(i)};    
end



