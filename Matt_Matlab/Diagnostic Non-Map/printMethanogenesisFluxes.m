function [rxns,fluxes] = printMethanogenesisFluxes(model,solution,map_flag)

%Take in the M. maripaludis model and a solution, then pull and print the
%fluxes of methanogenesis reactions
%Inputs
%model: the M. maripaludis model
%solution: a flux distribution solution to the M. maripaludis model
%map_flag: change to true (not a string)to draw the reactions on a Paint4Net map(default = false)

%Check for map flag
if (nargin<3)
    map_flag = false;
end

%Specify the reactions we want
rxns = {...
%Step 1: CO2 -> Formyl-MFR
'rxn11938_c0';...
%Step 2: Formyl-MFR -> Formyl-H4MPT
'rxn02431_c0';...
%Step 3: Formyl-H4MPT -> Methenyl-H4MPT
'rxn02480_c0';...
%Step 4: Methenyl-H4MPT -> Methylene-H4MPT
'rxn06696_c0';...
%Step 5: Methylene-H4MPT -> Methyl-H4MPT
'rxn03085_c0';...
%Step 6: Methyl-H4MPT -> Methyl-CoM
'rxn03020_c0';...
%Step 7: Methyl-CoM -> Methane
'rxn03127_c0';...
%Step 8: CoM-S-S-CoB -> CoM
'HdrABC';...
%Ferredoxin Generator
'Eha/Ehb';...
%ATP Synthase
'ATPS'};

%Create Reaction Strings to store reaction shorthand
rxn_strings = {...
'CO2 -> Formyl-MFR';...
'Formyl-MFR -> Formyl-H4MPT';...
'Formyl-H4MPT -> Methenyl-H4MPT';...
'Methenyl-H4MPT -> Methylene-H4MPT';...
'Methylene-H4MPT -> Methyl-H4MPT';...
'Methyl-H4MPT -> Methyl-CoM';...
'Methyl-CoM -> Methane';...
'CoM-S-S-CoB -> CoM';...
'Ferredoxin Generator';...
'ATP Synthase'};
%Intersect those reactions with the model, but keep it in the order of rxns
[~,~,idxB] = intersect(rxns,model.rxns,'stable');

%Use the index to pull out fluxes
fluxes = solution.x(idxB);

for i=1:length(rxns)
%Print it all out in pretty format
fprintf('\n%s(%s): %0.2f',rxn_strings{i},rxns{i},fluxes(i));
end
fprintf('\n\n')

if map_flag
    all_rxns = findMetSynPath(model,'EX_Methane_c0',solution,{'biomass0'},0.1);
    drawRxnsMap(model,all_rxns,solution);
end