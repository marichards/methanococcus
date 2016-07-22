function [id,name] = translateMetIDAndName(model,term)

%Given a metabolite ID or name, translate it to the complementary entity
%
% INPUT
% model: a COBRA Toolbox model structure
% term: either a metabolite ID or metabolite name from the supplied model
% to be translated
%
% OUTPUT
% id: the metabolite id for the supplied term
% name: the metabolite name for the supplied term
%
% Matthew Richards, 09/24/2015

% Check for whether it's an id
if ismember(term,model.mets)
    % Then assign the id and index
    [id,idx] = intersect(model.mets,term);
    % And use the index to get the name
    name = model.metNames{idx};

% If it's a name
elseif ismember(term,model.metNames)
    % Then assign the name and index
    [name,idx] = intersect(model.metNames,term);
    % And use the index to get the id
    id = model.mets{idx};           
end
        