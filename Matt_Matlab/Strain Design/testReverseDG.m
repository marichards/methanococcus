function gibbs = testReverseDG(model)
%%
% Takes in the M. maripaludis model and does a sensitivity analysis of
% concentration ratios for the main metabolites (CH4, CH3OH, CO2, H2O, H2,
% NH3, PO4)

% INPUT
% model: the M. maripaludis model, a COBRA toolbox model structure
%
% OUTPUT
%
% Matthew Richards, 10/13/2015
%%

% Set the substrate reactions
substrate_rxns = {...
    % Methane
    'EX_cpd01024[e0]';...
    % CO2
    'EX_cpd00011[e0]';...
    % Methanol
    'EX_cpd00116[e0]';...
    % H2 is not significant
    %'EX_cpd11640[e0]';...
    % NH3
    'EX_cpd00013[e0]';...
    % PO4 is not significant
    %'EX_cpd00009[e0]'...
    };

% Set the initial concentration vector
concentrations = ones(size(substrate_rxns));

% Do the initial calculation and grab the model with the correct bounds
[~,gibbs_flux(1),model] = reverseMethanogenesis(model,substrate_rxns,concentrations);

% Now the model has the dG reaction and proper fluxes, just do
% optimizeThermoModel() from this point forward

% Generate a set of log-spaced points for each reaction concentration
% Unlikely to go more than 10 M (10^4 mM) or less 10 nM (10-5 mM)
concs = logspace(-5,4,10);

% Initiate a matrix to hold gibbs fluxes
gibbs = zeros(length(concs),length(substrate_rxns));

% Grab the size of the gibbs matrix 
[magnitudes,rxns] = size(gibbs);

% Cycle through all the combinations and generate predictions of gibbs for
% each
for i = 1:rxns
    for j = 1:magnitudes
        
        % Alter the correct index 
        concentrations(i) = concs(j);
        
        % Simulate growth and pull out gibbs flux, store in a matrix
        [~,gibbs(j,i)] = optimizeThermoModel(model,substrate_rxns...
            ,concentrations,310,'EX_cpd00001[e0]');
        
        % Reset the concentrations
        concentrations = ones(size(substrate_rxns));
    end
end

% Plot a figure for each. Do subplots
rows = ceil(length(substrate_rxns)/2);
figure(1);
hold on;
for i = 1:length(substrate_rxns)
   
    subplot(rows,2,i)
    plot(concs,gibbs(:,i))
end




