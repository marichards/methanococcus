function omix_strings = formatRxnsForOmix(model,rxns,filename)

%For reactions, want to print it as such:
%reaction name : reaction formula 

%Make a filename
file_id = fopen(sprintf('%s.txt',filename),'w');

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
    %Print to a file
    fprintf(file_id,'%s: %s\n',model.rxnNames{rxn_idx},formula{1});
end