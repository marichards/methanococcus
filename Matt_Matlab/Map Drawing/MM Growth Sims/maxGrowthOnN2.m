function [solution,gibbs_flux,model] = maxGrowthOnN2(model,substrate_rxns,concentrations,print_flag)

% Simulate M. maripaludis growth on CO2 and H2 media, with nitrogen gas as
% the nitrogen source. Print out the growth rate and relevant fluxes,
% return the full solution, the predicted free energy generation, and the
% modified model with the overall Gibbs free energy reaction added to the S
% matrix
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
% 
% OPTIONAL INPUT
% substrate_rxns: a list of exchange reactions in the M. maripaludis model
% for which a known concentration will be supplied. If supplied, it must be
% accompanied by a corresponding "concentrations" array. (Default =
% {'EX_cpd00011[e]','EX_cpd11640[e0]','EX_cpd01024[e0]'})
% concentrations: a list of effective concentrations in mM corresponding to
% the exchange reactions listed in "substrate_rxns". (Default = [1 1 1])
%
% OUTPUT
% solution: a flux distribution solution from running FBA on the M.
% maripaludis model that maximizes biomass yield
% gibbs_flux: model prediction of overall free energy generation, based on
% the model exchange fluxes in the solution
% model: the M. maripaludis model, with an additional reaction
% (GIBBS_kJ/GDW) that predicts overall free energy generation
% 
% Matthew Richards, 09/24/2015

% Check if print_flag is supplied
if nargin < 4
    % Set default to true
    print_flag = true;
end

% Default model is growing on H2, but be sure it's only on H2.
model = switchToH2(model);

% Make sure that nitrogen is on, not ammonia or alanine
model = switchToN2(model);

% Specify substrate reactions and concentrations as 1 mM if not given
if nargin<2
    
    substrate_rxns = {'EX_cpd00011[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]'};
    concentrations = [1 1 1];
    warning_flag = 1;
end

%Solve by maximizing biomass
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]');
%Pull out the overall reaction CO2 + 4H2 --> CH4 + 2H2O

% Check for print flag 
if print_flag

    %Find the reaction indices
    [~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
    [~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
    [~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
    [~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
    [~,form_idx] = intersect(model.rxns,'EX_cpd00047[e0]');
    [~,n2_idx] = intersect(model.rxns,'EX_cpd00528[e0]');
    [~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');
    [~,ac_idx] = intersect(model.rxns,'EX_cpd00029[e0]');

    %Print the biomass flux
    fprintf('\n\nBiomass flux: %f\n\n',solution.f);
    %Print the reaction fluxes
    fprintf('Formate flux: %f\n',solution.x(form_idx))
    fprintf('CO2 flux: %f\n',solution.x(co2_idx))
    fprintf('H2 flux: %f\n',solution.x(h2_idx))
    fprintf('H2O flux: %f\n',solution.x(h2o_idx))
    fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
    fprintf('N2 flux: %f\n',solution.x(n2_idx))
    fprintf('PO4 flux: %f\n',solution.x(po4_idx))
    fprintf('Acetate flux: %f\n',solution.x(ac_idx))

    %Print the per-CO2 actual reaction
    fprintf('\nOverall reaction:\nCO2 + 4 H2 --> 2 H2O + CH4\n')
    fprintf('\nModel overall reaction (per mole CH4)\n')
    fprintf('%0.2f CO2 + %0.2f H2 --> %0.2f H2O + CH4\n\n',-solution.x(co2_idx)/solution.x(ch4_idx),...
        -solution.x(h2_idx)/solution.x(ch4_idx),solution.x(h2o_idx)/solution.x(ch4_idx))

    %Print the yield coefficient (grams biomass per mole CH4 produced)
    fprintf('Predicted Yield Coefficient: %0.2f gDCW/mol CH4\n\n',solution.f*1000/solution.x(ch4_idx)/log(2))

    %Find the ATP reaction index
    [~,atp_idx] = intersect(model.rxns,'ATPS');
    %Print the ATP yield coefficient (ATP per CH4)
    fprintf('Expected ATP/CH4 Yield: 0.5\n')
    fprintf('Predicted ATP/CH4 Yield: %0.3f\n\n', solution.x(atp_idx)/solution.x(ch4_idx))
end

%Add a warning for simulations with no concentrations given
if warning_flag
    warning('All external metabolite concentrations set to 1 mM');
end

% Print out the gibbs free energy prediction
if print_flag
    
    fprintf('Predicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)

end
