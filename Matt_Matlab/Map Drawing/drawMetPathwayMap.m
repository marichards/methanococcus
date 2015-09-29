function [rxns,rxnNames,fluxes] = drawMetPathwayMap(model,met,solution,biomass_rxn)

% Take in model, an FBA solution, and metabolite of interest. Pull the
% synthesis pathway for the metabolite, then draw the flux map fo the
% pathway using the Paint4Net toolbox.
%
% INPUT
% model:  a COBRA Toolbox model structure
% met: a metabolite of interest in the supplied model
% solution: a flux distribution solution to the supplied model
%
% OPTIONAL INPUT
% biomass_rxn: the name of the biomass reaction in the supplied model,
% which will be excluded from the metabolite synthesis pathway to greatly
% simplify the results. (Default = 'biomass0')
%
% OUTPUT
% rxns: a list of the reactions that are part of the synthesis pathway for
% the supplied metabolite of interest
% rxnNames: a list of reaction names for the reactions in "rxns"
% fluxes: a list of fluxes in the supplied solution corresponding to the
% reactions in "rxns"
% 
% Matthew Richards, 09/28/2015


% Check for a biomass reaction
if (nargin<4)
    biomass_rxn = 'biomass0';
end
% Use the findMetSynPath code with its default values to find the
% metabolite synthesis pathway, excluding the biomass reaction
[rxns,rxnNames,fluxes] = findMetSynPath(model,met,solution,{biomass_rxn});

% Use Paint4Net to draw the reaction map as a biograph
draw_by_rxn(model,rxns,'true','struc',{''},{''},solution.x)
