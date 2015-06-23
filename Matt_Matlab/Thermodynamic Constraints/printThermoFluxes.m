function printThermoFluxes(model,solution)

%Written by Matthew Richards 03/01/2013
%
%Input: model with overallDG and freeEnergy parameters, solution
%
%Output: exchanges, including overallDG and biomass

%Find them (non-zero in freeEnergy)
exc_rxns = find(model.freeEnergy~=0);
%Add biomass and overall DG


labels = model.rxns(exc_rxns);
fluxData = solution.x(exc_rxns,:);
printLabeledData(labels,fluxData)
    
    