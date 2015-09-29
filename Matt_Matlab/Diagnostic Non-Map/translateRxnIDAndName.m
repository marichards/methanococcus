function [id,name] = translateRxnIDAndName(model,term)

%Given a reaction ID or name, translate it to the complementary entity
%
% INPUT
% model: a COBRA Toolbox model structure
% term: either a reaction ID or reaction name from the supplied model to be
% translated 
%
% OUTPUT
% id: the reaction id for the supplied term
% name: the reaction name for the supplied term
%
% Matthew Richards, 09/24/2015

% Check for whether it's an id
if ismember(term,model.rxns)
    % Then assign the id and index
    [id,idx] = intersect(model.rxns,term);
    % And use the index to get the name
    name = model.rxnNames{idx};

% If it's a name
elseif ismember(term,model.rxnNames)
    % Then assign the name and index
    [name,idx] = intersect(model.rxnNames,term);
    % And use the index to get the id
    id = model.rxns{idx};    
end
        