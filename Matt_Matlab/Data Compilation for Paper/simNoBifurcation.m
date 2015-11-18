function solution = simNoBifurcation(model)

% Simulate M. maripaludis growth with ferredoxin reduction removed from HDR. Do so on CO2 and H2 media, with ammonia as the
% nitrogen source. Print out the growth rate and relevant fluxes, return
% the full solution.
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
% 
% OUTPUT
% solution: a flux distribution solution from running FBA on the M.
% maripaludis model that maximizes biomass yield
% 
% Matthew Richards, 11/17/2015

% Modify the HDR reaction to not include ferredoxin:

model = addReaction(model,{'HdrABC','HdrABC'},...
    'cpd02935[c0] + cpd11640[c0] -> cpd02246[c0] + cpd02817[c0]');

% Simulate model growth on H2
solution = maxGrowthOnH2(model);