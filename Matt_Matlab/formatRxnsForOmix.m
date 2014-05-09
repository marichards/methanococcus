function omix_strings = formatRxnsForOmix(model,rxns)

%For reactions, want to print it as such:
%reaction name : reaction formula 

%Create an array for strings
omix_strings = cell(length(rxns),1);

for i=1:length(rxns)
%Grab the formula
formula = printRxnFormula(model,rxns{i},'False','False');
%Grab the name of the reaction
%First grab the index
[~,rxn_idx] = intersect(model.rxns,rxns{i});
%Now use the index to grab the name and make a string
omix_strings{i}=sprintf('%s: %s',model.rxnNames{rxn_idx},formula{1});
end