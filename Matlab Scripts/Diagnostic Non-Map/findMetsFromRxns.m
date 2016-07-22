function mets = findMetsFromRxns(model,rxns)

% A parallel function to "findRxnsFromMets" that takes in a list of
% reactions and returns a list of all metabolites in the given reactions
%
% INPUT:
% model: a COBRA Toolbox model structure
% rxns: a list of reactions to look for in the supplied model
%
% OUTPUT:
% mets: the list of metabolites that appear in any of the 'rxns'
%
% Matthew Richards, 09/24/2015


% First find the reaction index
[~,idx] = intersect(model.rxns,rxns);

%Using that index, find all indices of metabolites (less than 1)
mets= {};
for i=1:length(idx)
mets = [mets; model.mets(model.S(:,idx(i))~=0)];
end
%Trim non-unique metabolites
mets=unique(mets);
end