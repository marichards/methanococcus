function createGeneRxnDict(model,workbook,sheet)

%This code will pull out all reactions in the model and their gene
%associations, then write those to an excel sheet
%
%Inputs:
%model: a cobra model structure
%workbook: an excel workbook to write to (.xlsx)
%sheet: an excel sheet to write to

formulas = printRxnFormula(model,model.rxns,'False');
A = [model.rxns,model.rxnNames,model.grRules,formulas,num2cell(model.lb),num2cell(model.ub)];
xlswrite(workbook,A,sheet);