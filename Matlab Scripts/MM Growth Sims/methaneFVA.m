function [fva_min,fva_max] = methaneFVA(model,optPercent)

% Takes in the M. maripaludis model, runs flux variability analysis to the
% specified percentage, and then prints out the major chemical fluxes. It
% also returns the minimum and maximum fluxes for all reactions. 
%
%%%% NOTE: There's a problem with the general fluxVariability code in the
%%%% COBRA toolbox, so if you get an error about "Index exceeds matrix
%%%% dimensions", you should first try just re-running it.
%
% INPUT 
% model: the M. maripaludis model, a COBRA toolbox model structure
% 
% OPTIONAL INPUT
% optPercent: a value from 0-100 that specifies by the minimum acceptable
% value of the objective function by percentage of the maximum value
% achievable. (Default = 100)
% 
% OUTPUT
% fva_min: an array of the minimum flux bound corresponding to every
% reaction in the M. maripaludis model
% fva_max: an array of the maximum flux bound corresponding to every
% reaction in the M. maripaludis model
%
% Matthew Richards, 09/24/2015

% If no optPercent is specified, set it to 100%
if (nargin <2)
    optPercent = 100;
end

%First run FVA
[fva_min,fva_max] = fluxVariability(model,optPercent);

%Pull out the indices of important things:
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
[~,form_idx] = intersect(model.rxns,'EX_cpd00047[e0]');
[~,bio_idx] = intersect(model.rxns,'biomass0');

%Print out the flux range for each one 
fprintf('\n\nBiomass flux range: %0.4f-%0.4f\n',fva_min(bio_idx),fva_max(bio_idx));
fprintf('CO2 flux: %0.2f-%0.2f\n',fva_min(co2_idx),fva_max(co2_idx))
fprintf('H2 flux: %0.2f-%0.2f\n',fva_min(h2_idx),fva_max(h2_idx))
fprintf('H2O flux: %0.2f-%0.2f\n',fva_min(h2o_idx),fva_max(h2o_idx))
fprintf('CH4 flux: %0.2f-%0.2f\n',fva_min(ch4_idx),fva_max(ch4_idx))
fprintf('Formate flux range: %0.2f-%0.2f\n\n',fva_min(form_idx),fva_max(form_idx))