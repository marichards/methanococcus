function printRxnFluxes(model,solution,rxns)

% Take in a model and a flux solution, plus a list of reactions we're
% interested in. Go into the flux solution and print out the flux going
% through the specified reactions
%
% INPUT
% model: a COBRA Toolbox model structure
% solution: a flux distribution solution to the supplied model
% rxns: a list of reactions of interest in the supplied model
%
% Matthew Richards, 09/24/2015


%Find the indices of the reactions
[rxns,idx] = intersect(model.rxns,rxns);

%Pull out the proper fluxes
fluxes = solution.x(idx);

%Print out the fluxes for all the reactions
fprintf('\n')
for i=1:length(rxns)
    fprintf('%s\t%0.3f\n',rxns{i},fluxes(i))
end
