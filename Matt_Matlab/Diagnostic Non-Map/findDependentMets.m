function mets = findDependentMets(model,rxn)

% Pulls out a reaction, tests essentiality of the reaction, and if it's
% essential then finds what mets depend on its inclusion

% INPUT
% model: a COBRA Toolbox model structure
% rxn: a reaction in the supplied model to evaluate
%
% OUTPUT
% mets: a list of biomass metabolites that depend on the reaction
%
% Matthew Richards, 09/24/2015


%Test essentiality of the reaction by removing it and simulating growth
model = removeRxns(model,rxn);
solution = optimizeCbModel(model);

%If the solution.f is effectively 0, then do the biomass precursor check
if solution.f < 1e-10
    mets = biomassPrecursorCheck(model);
else
    mets = {};
end

    