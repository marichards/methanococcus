function rxns = printRxnsFromGenes(model,genes)

%rxns = printRxnsFromGenes(model,genes)
%
%Inputs
%
%Outputs

%Use the built-in function to grab the reactions
rxn_dict = findRxnsFromGenes(model,genes);

%Grab and cycle through the structure fields to make an array of reactions
fields = fieldnames(rxn_dict);

rxns={};
for i=1:numel(fields)
m=size(rxn_dict.(fields{i}));
total = total+m;
end

%Now do it again but save the reactions
