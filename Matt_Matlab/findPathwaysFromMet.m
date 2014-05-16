function findPathwaysFromMet(model,solution,met_name)

%Given a metabolite name, a model, and a solution, find the reactions that
%contain the metabolite AND have non-zero fluxes.  Print out the reaction
%names with their formulas, then print out the fluxes

%First use pull_rxns to get the name of the reactions and print the
%formulas
rxns = findRxnsFromMets(model,met_name);

%Now find the indices for those reactions
[rxns,idx]=intersect(model.rxns,rxns);

%Find the reactions in the solution and, if they have non-zero fluxes,
%print them out with the fluxes
for i=1:length(idx)
   
    if solution.x(idx(i))~= 0 
        
        fprintf('%s flux: %f\n',rxns{i},solution.x(idx(i)));
    end
    
end