function printRxnFluxes(model,solution,rxns)

%Take in a model and a flux solution, plus a list of reactions we're
%interested in. Go into the flux solution and pull out the flux going
%through the specified reactions

%Find the indices of the reactions
[rxns,idx] = intersect(model.rxns,rxns);

%Pull out the proper fluxes
fluxes = solution.x(idx);

%Print out the fluxes for all the reactions
fprintf('\n')
for i=1:length(rxns)
    fprintf('%s\t%0.3f\n',rxns{i},fluxes(i))
end
