function [solution,gibbs_flux,model] = maxGrowthOnMethanolOnly(model,substrate_rxns,concentrations)

% Simulate a growth of a mutated M. maripaludis strain with the ability to
% uptake and utilize Methanol + H2 for growth and methanogenesis, with
% ammonia as the nitrogen source. Print out the growth rate and relevant
% fluxes, return the full solution, the predicted free energy generation,
% and the modified model with the overall Gibbs free energy reaction added
% to the S matrix
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
% 
% OPTIONAL INPUT
% substrate_rxns: a list of exchange reactions in the M. maripaludis model
% for which a known concentration will be supplied. If supplied, it must be
% accompanied by a corresponding "concentrations" array. (Default =
% {'EX_cpd00116[e]','EX_cpd11640[e0]','EX_cpd01024[e0]'})
% concentrations: a list of effective concentrations in mM corresponding to
% the exchange reactions listed in "substrate_rxns". (Default = [1 1 1])
%
% OUTPUT
% solution: a flux distribution solution from running FBA on the M.
% maripaludis model that maximizes biomass yield
% gibbs_flux: model prediction of overall free energy generation, based on
% the model exchange fluxes in the solution
% model: the M. maripaludis model, with the added methanol uptake pathway
% and an additional reaction (GIBBS_kJ/GDW) that predicts overall free
% energy generation. Note that because the first step of hydrogenotrophic
% methanogenesis has been constrained to 0 in this model, it will not grow
% hydrogenotrophically unless this constraint is removed. 
% 
% Matthew Richards, 09/28/2015

% Add in the reaction for converting methanol to methyl-coM
model = addReaction(model,{'rxn10568[c0]','Methanol: coenzyme M methyltransferase'},...
    'cpd00116[c0] + cpd02246[c0] <=> cpd02438[c0] + cpd00001[c0]');

% Add in the methanol uptake
model = addReaction(model,{'rxn10570[c0]','Methanol diffusion'},...
    'cpd00116[e0] <=> cpd00116[c0]');
model = addReaction(model,{'EX_cpd00116[e0]','EX_Methanol[e0]'},...
    'cpd00116[e0] <=> ');

% Add metabolite info for methanol
[~,idx] = intersect(model.mets,'cpd00116[c0]');
model.metNames{idx} = 'Methanol[c0]';
model.metCharge(idx)=0;
model.metFormulas{idx}='CH4O';
[~,idx] = intersect(model.mets,'cpd00116[e0]');
model.metNames{idx} = 'Methanol[e0]';
model.metCharge(idx)=0;
model.metFormulas{idx}='CH4O';

% Add it to the free energy
[~,idx] = intersect(model.rxns,'EX_cpd00116[e0]');
model.freeEnergy(idx) = -0.0302;

% Make certain that H2 is on and acetate/formate are off
% Turn H2 off
model = changeRxnBounds(model,'EX_cpd11640[e0]',0,'l');
% Turn down formate
model = changeRxnBounds(model,'EX_cpd00047[e0]',0,'l');
% Turn down acetate
model = changeRxnBounds(model,'EX_cpd00029[e0]',0,'l');
% Make certain that NH3 is on, N2 is off, Alanine is off
% Turn up ammonia
model = changeRxnBounds(model,'EX_cpd00013[e0]',-1000,'l');
% Turn down alanine
model = changeRxnBounds(model,'EX_cpd00035[e0]',0,'l');
% Turn down nitrogen
model = changeRxnBounds(model,'EX_cpd00528[e0]',0,'l');

% Turn off CO2
model = changeRxnBounds(model,'EX_cpd00011[e0]',0,'l');

% Set a bound on methane
model = changeRxnBounds(model,'EX_cpd01024[e0]',46,'b');

% Check if model is specific ferredoxins or not and set bound on Eha and
% Ehb in either case
if ismember('Eha/Ehb',model.rxns)
    % If not, then set bounds on Eha/Ehb
    model = changeRxnBounds(model,'Eha/Ehb',46,'u');
    model = changeRxnBounds(model,'Eha/Ehb',-46,'l');
else
    % If it is, then sent on both Eha and Ehb
    model = changeRxnBounds(model,'Eha',46,'u');
    model = changeRxnBounds(model,'Eha',-46,'l');
    model = changeRxnBounds(model,'Ehb',46,'u');
    model = changeRxnBounds(model,'Ehb',-46,'l');
end

% Specify substrate reactions and concentrations as 1 mM if not given
if nargin<2    
    substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]'};
    concentrations = [1 1 1];
    warning_flag = 1;
end
% Solve by maximizing biomass
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]');

% Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(model.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');


% Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);

if solution.f > 0 
% Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))

% Print the per-CO2 actual reaction
fprintf('\nOverall reaction:\nCH3OH + H2 --> H2O + CH4\n')
fprintf('\nModel overall reaction (per mole CH4)\n')
fprintf('%0.2f CH3OH + %0.2f H2 --> %0.2f H2O + CH4\n\n',-solution.x(meoh_idx)/solution.x(ch4_idx),...
    -solution.x(h2_idx)/solution.x(ch4_idx),solution.x(h2o_idx)/solution.x(ch4_idx))

% Print the yield coefficient (grams biomass per mole CH4 produced)
fprintf('Predicted Yield Coefficient: %0.2f gDCW/mol CH4\n\n',solution.f*1000/solution.x(ch4_idx)/log(2))

% Find the ATP reaction index
[~,atp_idx] = intersect(model.rxns,'ATPS');
% Print the ATP yield coefficient (ATP per CH4)
fprintf('Expected ATP Yield: 0.5\n')
fprintf('Predicted ATP Yield: %0.3f\n\n', solution.x(atp_idx)/solution.x(ch4_idx))

% Print out the gibbs free energy prediction
% Add a warning for simulations with no concentrations given
if warning_flag
    warning('All external metabolite concentrations set to 1 mM');
end
fprintf('Predicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)

end
end