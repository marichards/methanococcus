function [gam,ngam,atp_fluxes] = determineATPM(model,growth_rates,ch4_rates)

% Using measured experimental growth rates and methane evolution rate,
% constrain the model to measured rates and measure ATP hydrolysis that
% results. Repeat and graph ATPM vs growth rate to calculcate ATP
% maintenance

% Check to make sure growth rates and ch4 rates are the same length
if length(growth_rates) ~= length(ch4_rates)
    error('Please ensure that there are an equal number of growth rates and methane evolution rates');
end

% Make an ATP-hydrolyzing model by removing ATP from biomass, adjusting
% bounds on NGAM, and making NGAM the objective function

% First, find ATP
[~,idx] = intersect(model.mets,'cpd00002[c0]');
% Find index of biomass
[~,bio_idx] = intersect(model.rxns,'biomass0');
% Remove it from biomass
model.S(idx,bio_idx) = 0;
% Repeat for ADP and Phosphate
[~,idx] = intersect(model.mets,'cpd00008[c0]');
% Find index of biomass
[~,bio_idx] = intersect(model.rxns,'biomass0');
% Remove it from biomass
model.S(idx,bio_idx) = 0;
[~,idx] = intersect(model.mets,'cpd00009[c0]');
% Find index of biomass
[~,bio_idx] = intersect(model.rxns,'biomass0');
% Remove it from biomass
model.S(idx,bio_idx) = 0;

% Change the objective to the NGAM reaction and change its bounds
model = changeObjective(model,'biomass0',0);
model = changeObjective(model,'rxn00062[c0]',1);
model = changeRxnBounds(model,'rxn00062[c0]',inf,'u');
model = changeRxnBounds(model,'rxn00062[c0]',-inf,'l');

% Initiate the atp_flux vector
atp_fluxes = zeros(size(growth_rates));

% Iterate through the supplied list of growth rates and ch4_rates 
for i = 1:length(growth_rates)
    
    % Constrain the model by the supplied growth rate and ch4 rate
    model = changeRxnBounds(model,'biomass0',growth_rates(i),'b');
    model = setMethaneSecretion(model,ch4_rates(i));    
    
    % Simulate the model
    solution = optimizeCbModel(model,[],'one');

    % Pull out the flux of ATP and add to vector of ATP fluxes
    atp_fluxes(i) = solution.f;
end

% Pull out the slope and intercept of the atp fluxes vs growth rates
p = polyfit(growth_rates,atp_fluxes,1);

gam = p(1);
ngam = p(2);


