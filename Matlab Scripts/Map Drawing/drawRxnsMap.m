function drawRxnsMap(model,rxns,solution)

% A small shortcut to use the Paint4Net toolbox. Takes in a model, a set of
% reactions of interest, and an FBA solution. Uses the Paint4Net function
% draw_by_rxn to create a flux map of the supplied reactions. 
%
% INPUT
% model:  a COBRA Toolbox model structure
% met: a list of reactions of interest in the supplied model
% solution: a flux distribution solution to the supplied model
% 
% Matthew Richards, 09/28/2015


%Take in a model, a reaction set of interest, and flux solution, then draw a map with Paint4Net
draw_by_rxn(model,rxns,'true','struc',{''},{''},solution.x);
