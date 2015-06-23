function printReverseFluxes(model,solution,map_flag)

% For the reverse methanogenesis maripaludis model, print out the main
% fluxes so I can see what's happening

%Check for map flag
if (nargin<3)
    map_flag = false;
end

%Specify the reactions we want
rxns = {...
%Step 1: Methane -> Methyl-CoM
'rxn03127_c0';...
%Step 2: Methyl-CoM -> Methanol
'rxn10568_c0';...
%Step 3: CoM-S-S-CoB -> CoM
'HdrABC';...
%Step 4: Na+ Generator
'Eha/Ehb';...
%ATP Synthase
'ATPS'};

%Create Reaction Strings to store reaction shorthand
rxn_strings = {...
'Methane -> Methyl-CoM';...
'Methyl-CoM -> Methanol';...
'CoM-S-S-CoB -> CoM';...
'Ferredoxin Generator';...
'ATP Synthase'};
%Intersect those reactions with the model, but keep it in the order of rxns
[~,~,idxB] = intersect(rxns,model.rxns,'stable');

%Print out the reaction formulas
printRxnFormula(model,rxns);

%Use the index to pull out fluxes
fluxes = solution.x(idxB);

for i=1:length(rxns)
%Print it all out in pretty format
fprintf('\n%s(%s): %0.2f',rxn_strings{i},rxns{i},fluxes(i));
end
fprintf('\n\n')

if map_flag
    all_rxns = findMetSynPath(model,'Methanol_c0',solution,{'biomass0'},0.1);
    drawRxnsMap(model,all_rxns,solution);
end