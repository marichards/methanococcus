function mets = findMetsFromRxns(model,rxns)

%Find the metabolites in a given reaction

%Inputs:
%model: a COBRA model structure
%rxns: a list of reactions to look for

%Outputs:
%mets: the list of metabolites that appear in any of the 'rxns'

%First find the reaction index
[~,idx] = intersect(model.rxns,rxns);

%Using that index, find all indices of metabolites (less than 1)
mets= {};
for i=1:length(idx)
mets = [mets; model.mets(model.S(:,idx(i))~=0)];
end
%Trim non-unique metabolites
mets=unique(mets);
end