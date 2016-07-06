function demoFreeEnergy(model)

% Demonstrate the capabilities of calculating overall free energy by
% varying hydrogen concentrations and plotting those concentrations against
% their associated gibbs free energies
%
% INPUT
% model: the M. maripaludis model, a COBRA toolbox model
%
% OUTPUT

% Set x on a log scale from 10^-10 to 1
h2_concs = logspace(-10,0);

% Create an array of gibbs free energies corresponding to the H2
% concentrations
gibbs = zeros(length(h2_concs),1);
for i = 1:length(h2_concs)
[~,gibbs(i)] = maxGrowthOnH2Only(model,{'EX_cpd11640[e0]'},h2_concs(i),false);
end

% Plot the results on a semilog plot
figure(1)
h = semilogx(h2_concs,gibbs);
set(h,'LineWidth',3)
xlabel('H_{2} Concentration (mM)','FontSize',14)
ylabel('Predicted \DeltaG (kJ/GDW)','FontSize',14)