function [rxns,fluxes,rxnNames] = drawRxnsMap(model,rxns,solution)

%Take in a model, a reaction set of interest, and flux solution, then draw a map with Paint4Net
draw_by_rxn(model,rxns,'true','struc',{''},{''},solution.x);
