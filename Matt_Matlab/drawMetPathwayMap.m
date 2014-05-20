function [rxns,fluxes,rxnNames] = drawMetPathwayMap(model,met,solution)

%Take in a model, a metabolite of interest, and flux solution, then draw a map with Paint4Net
[rxns,fluxes,rxnNames] = traceFullMetPath(model,met,solution,{'biomass0'});
draw_by_rxn(model,rxns,'true','struc',{''},{''},solution.x)
