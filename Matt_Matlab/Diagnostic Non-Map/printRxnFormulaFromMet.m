function printRxnFormulaFromMet(model,met)

% Take in a metabolite from a COBRA model, then print reaction formulas for
% all reactions with the metabolite

printRxnFormula(model,findRxnsFromMets(model,met));