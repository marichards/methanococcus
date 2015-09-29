function printRxnFormulaFromMet(model,met)

% A little shortcut that takes in a metabolite from a COBRA model, then
% print reaction formulas for all reactions that contain the metabolite
%
% INPUT
% model: a COBRA Toolbox model structure
% met: a metabolite of interest from the supplied model
%
% Matthew Richards, 09/24/2015

printRxnFormula(model,findRxnsFromMets(model,met));