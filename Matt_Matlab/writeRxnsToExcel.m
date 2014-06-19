function writeRxnsToExcel(model,rxns,workbook,sheet)

%This code will simulate a COBRA model for maximizing biomass, then write
%the reactions and fluxes to an Excel sheet. 

%Inputs:
%model: a cobra model structure
%workbook: an excel workbook to write to (.xlsx)
%sheet: an excel sheet to write to

%Solve the model, write the reactions, name, formula to an excel file
[rxns,idx]=intersect(model.rxns,rxns);
names = model.rxnNames(idx);
formulas = printRxnFormula(model,rxns,'False');
A = [rxns,names,formulas];
xlswrite(workbook,A,sheet);