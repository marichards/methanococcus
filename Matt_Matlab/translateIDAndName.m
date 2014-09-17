function [id,name] = translateIDAndName(model,term)

%Given a reaction ID or name, translate it to the other one
%
%Input
%model: a COBRA model structure
%term: either a reaction ID or name that you want translated
%
%Output
%id: the reaction id
%name: the reaction name

%Check for whether it's an id
if ismember(term,model.rxns)
    %Then assign the id and index
    [id,idx] = intersect(model.rxns,term);
    %And use the index to get the name
    name = model.rxnNames{idx};

%If it's a name
elseif ismember(term,model.rxnNames)
    %Then assign the name and index
    [name,idx] = intersect(model.rxnNames,term);
    %And use the index to get the id
    id = model.rxns{idx};       
    
end
        