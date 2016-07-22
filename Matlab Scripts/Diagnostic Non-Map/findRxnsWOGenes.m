function rxns = findRxnsWOGenes(model)

% Go through a model and find all reactions without genes (the gapfilled
% reactions). Don't return the exchange reactions, as those are model
% artifacts and shouldn't be considered the same as internal reactions
% without genes
%
% INPUT
% model: a COBRA Toolbox model structure
%
% OUTPUT
% rxns: the set of reactions from the model that lack genes, corresponding
% to the "gapfilled" reactions
%
% Matthew Richards, 09/24/2015

rxns = {};

for i =1:length(model.rxns)   
    % If there are no genes
    if isempty(model.grRules{i})        
        % Then save them
        rxns = [rxns;model.rxns{i}]; 
    end    
end

% Pull out exchanges
rxns = setdiff(rxns,model.rxns(findExcRxns(model)));