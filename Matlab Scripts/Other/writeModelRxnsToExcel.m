function writeModelRxnsToExcel(model,workbook,sheet)

%This code will simulate a COBRA model for maximizing biomass, then write
%the reactions and fluxes to an Excel sheet. 
%%% This is obsolete; just use "writeRxnsToExcel" with "model.rxns" as
%%% input

%Inputs:
%model: a cobra model structure
%workbook: an excel workbook to write to (.xlsx)
%sheet: an excel sheet to write to

%Solve the model, write the reactions, fluxes, ub/lb to an excel file

formulas = printRxnFormula(model,model.rxns,false);
A = [model.rxns,model.rxnNames,model.grRules,formulas,num2cell(model.lb),num2cell(model.ub)];
xlswrite(workbook,A,sheet);