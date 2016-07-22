function essentials = findEssentialRxns(model,rxns,cutoff)

% Based on a cutoff for percentage of growth considered "lethal", find the
% essential reactions in the given list of reactions
%
% INPUT
% model: a COBRA Toolbox model structure
%
% OPTIONAL INPUT
% rxns: a list of reactions in the supplied model to evaluate for 
% essentiality (Default = model.rxns)
% cutoff: a supplied cutoff for the ratio of knockout growth to wild-type
% growth that signifies lethality. Scales from 0-1. (Default = 0)
% 
% OUTPUT
% essentials: subset of supplied "rxns" that are essential for growth
% based upon the supplied cutoff
%
% Matthew Richards, 09/24/2014

%Qualify if there isn't a cutoff; make it 0
if (nargin<3)
    cutoff = 0;
end

%If there's no genes given, then just use all of them
if (nargin<2)
    rxns = model.rxns;
end

% Run the single gene deletion
grRatio = singleRxnDeletion(model,'FBA',rxns);

%Pull out things where the growth ratio is below the cutoff
essentials = rxns(grRatio<=cutoff);