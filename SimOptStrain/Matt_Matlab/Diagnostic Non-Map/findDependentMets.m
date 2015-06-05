function mets = findDependentMets(model,rxn)

%Pull out a reaction, test essentiality of the reaction, and if it's
%essential then find what mets depend on its inclusion

%Inputs
%model: a COBRA model structure
%rxn: a reaction to evaluate

%Outputs
%mets: a list of biomass metabolites that depend on the reaction

%%%%%%%%%%%%%%%

%Test essentiality of the reaction by removing it and simulating growth
model = removeRxns(model,rxn);
solution = optimizeCbModel(model);

%If the solution.f is 0, then do the biomass precursor check
if solution.f == 0
    mets = biomassPrecursorCheck(model);
else
    mets = {};
end

    