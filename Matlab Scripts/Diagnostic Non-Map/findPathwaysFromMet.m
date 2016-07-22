function findPathwaysFromMet(model,solution,met)

% Given a metabolite name, a model, and a solution, find the reactions that
% contain the metabolite AND have non-zero fluxes.  Print out the reaction
% names with their formulas, then print out the fluxes
%
% INPUT
% model: a COBRA Toolbox model structure
% solution: a flux distribution solution for the supplied model
% met_name: a metabolite of interest in the supplied model
%
% Matthew Richards, 09/24/2015


% First find all reactions that contain the metabolite of interest
rxns = findRxnsFromMets(model,met);

% Now find the indices for those reactions
[rxns,idx]=intersect(model.rxns,rxns);

% Find the reactions in the solution and, if they have non-zero fluxes,
% print them out with the fluxes
for i=1:length(idx)   
    if solution.x(idx(i))~= 0         
        fprintf('%s flux: %f\n',rxns{i},solution.x(idx(i)));
    end    
end