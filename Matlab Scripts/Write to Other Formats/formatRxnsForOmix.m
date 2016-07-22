function omix_strings = formatRxnsForOmix(model,rxns,filename)
%%
%Formats a supplied set of reactions into a text file for import into the
%network visualization software, OMIX.
%
% INPUT
% model: a COBRA Toolbox model structure
% rxns: a cell array containing reaction IDs of interestfor the supplied
% model
% filename: a character string used to name the output file
% 
% Matthew Richards, 09/29/2015


%%
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