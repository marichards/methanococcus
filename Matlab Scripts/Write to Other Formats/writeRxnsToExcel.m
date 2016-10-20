function writeRxnsToExcel(model,rxns,workbook,sheet)
%%
% Takes a specified set of reactions from a COBRA model and writes them to
% an Excel file, along with their names, GPR relationships, and
% formulas.
%
% INPUT:
% model: a COBRA Toolbox model structure
% rxns: a list of reaction IDs for reactions of interest
% workbook: an excel workbook to write to (.xlsx)
% sheet: a specified excel sheet to write to
%
% Matthew Richards, 09/29/2015


%%
% Write the reactions, name, genes, and formula to an excel file
[rxns,idx]=intersect(model.rxns,rxns);
names = model.rxnNames(idx);
genes = model.grRules(idx);
formulas = printRxnFormula(model,rxns,false,false,true);
A = [rxns,names,genes,formulas];
xlswrite(workbook,A,sheet);