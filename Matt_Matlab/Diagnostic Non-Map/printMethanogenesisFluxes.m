function [rxns,fluxes] = printMethanogenesisFluxes(model,solution,map_flag)

% Take in the M. maripaludis model and an FBA solution, then pull and print
% the fluxes of methanogenesis reactions. Optionally, also use Paint4Net to
% draw a flux map of the methanogenesis reactions.
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
% solution: a flux distribution solution to the M. maripaludis model
%
% OPTIONAL INPUT
% map_flag: a bianary variable that, when set to true, uses Paint4Net to
% draw a flux map of the methanogenesis reactions. (Default = false)
%
% OUTPUT
% rxns: a list of the reactions that are part of methanogenesis in M.
% maripaludis, aka "The Wolfe Cycle"
% fluxes: a list of fluxes corresponding to the methanogenesis reactions in
% "rxns"
%
% Matthew Richards, 09/24/2015


% Check for map flag
if (nargin<3)
    map_flag = false;
end

% Specify the reactions we want
rxns = {...
% Step 1: CO2 -> Formyl-MFR
'rxn11938[c0]';...
% Step 2: Formyl-MFR -> Formyl-H4MPT
'rxn02431[c0]';...
% Step 3: Formyl-H4MPT -> Methenyl-H4MPT
'rxn02480[c0]';...
% Step 4: Methenyl-H4MPT -> Methylene-H4MPT
'rxn06696[c0]';...
% Step 5: Methylene-H4MPT -> Methyl-H4MPT
'rxn03085[c0]';...
% Step 6: Methyl-H4MPT -> Methyl-CoM
'rxn03020[c0]';...
% Step 7: Methyl-CoM -> Methane
'rxn03127[c0]';...
% Step 8: CoM-S-S-CoB -> CoM
'HdrABC';...
% Ferredoxin Generator
'Eha/Ehb';...
% ATP Synthase
'ATPS'};

% Create Reaction Strings to store reaction shorthand
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
% Intersect those reactions with the model, but keep it in the order of rxns
[~,~,idxB] = intersect(rxns,model.rxns,'stable');

% Use the index to pull out fluxes
fluxes = solution.x(idxB);

for i=1:length(rxns)
% Print it all out in pretty format
fprintf('\n%s(%s): %0.2f',rxn_strings{i},rxns{i},fluxes(i));
end
fprintf('\n\n')

if map_flag
    all_rxns = findMetSynPath(model,'cpd01024[e0]',solution,{'biomass0'},0.1);
    drawRxnsMap(model,all_rxns,solution);
end