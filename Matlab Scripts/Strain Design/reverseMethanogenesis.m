function [solution,gibbs_flux,model] = reverseMethanogenesis(model,substrate_rxns,concentrations)

% Simulate a growth of a mutated M. maripaludis strain with the ability to
% uptake and utilize methane for growth and reverse methanogenesis to
% produce methanol with ammonia as the nitrogen source. Print out the
% growth rate and relevant fluxes, return the full solution, the predicted
% free energy generation, and the modified model with the overall Gibbs
% free energy reaction added to the S matrix
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
% solution: a flux distribution solution from running FBA on the modified
% M. maripaludis model that maximizes biomass yield
% gibbs_flux: model prediction of overall free energy generation, based on
% the model exchange fluxes in the solution
% model: the M. maripaludis model, with the added methanol uptake pathway,
% an added electron transfer from H2 to ferredoxin, and an additional
% reaction (GIBBS_kJ/GDW) that predicts overall free energy generation.
% Note that because the first step of hydrogenotrophic methanogenesis has
% been constrained to 0 in this model, it will not grow
% hydrogenotrophically unless this constraint is removed.
% 
% Matthew Richards, 09/28/2015


% First add in methanol pathway
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
[~,idx] = intersect(model.rxns,'rxn10570[c0]');
model.freeEnergy(idx) = 0;

% Change bounds to function in reverse
% Change bounds such that methane goes IN instead of OUT
model = changeRxnBounds(model,'EX_cpd01024[e0]',-1000,'l');
model = changeRxnBounds(model,'EX_cpd01024[e0]',0,'u');
% Force methanol to come OUT
model = changeRxnBounds(model,'EX_cpd00116[e0]',10,'l');
model = changeRxnBounds(model,'EX_cpd00116[e0]',10,'u');
% Turn off Hydrogen input and let it come out
model = changeRxnBounds(model,'EX_cpd11640[e0]',0,'l');
model = changeRxnBounds(model,'EX_cpd11640[e0]',1000,'u');
% Try turning off proton output, else it pours out. 
model = changeRxnBounds(model,'EX_cpd00067[e0]',0,'b');

% Be sure formate/acetate and N2/Alanine are off, and that NH3 is on
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

% Turn off FWD too, so we can't grow on CO2
model = changeRxnBounds(model,'rxn11938[c0]',0,'b');

% Add in the H2 -> Fd_red reactions 
model = addReaction(model,{'rxn05759[c0]','Reduced ferredoxin:H+ oxidoreductase'},...
    '2 cpd00067[c0] + cpd11620[c0] <=> cpd11640[c0] + cpd11621[c0]');
    
% Add a row to freeEnergy to make dimensions correct
[~,idx] = intersect(model.rxns,'rxn05759[c0]');
model.freeEnergy(idx) = 0;
% Specify substrate reactions and concentrations as 1 mM if not given
if nargin<2    
    substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]'};
    concentrations = [1 1 1];
    warning_flag = 1;
else
    warning_flag = 0;
end

% Remove the bounds on Eha
model = removeEhaBounds(model);

% Solve by maximizing biomass
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',false);

%Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(model.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');

if solution.f > 0 
%Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);
%Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))

%Print the per-CO2 actual reaction
fprintf('\nExpected Overall reactions:\nCH4 + H2O -> CH3OH + H2\n3 CH4 + CO2 + 2 H2O --> 4 CH3OH\n')
fprintf('\nModel overall reaction (per mole CH4):\n')
fprintf('CH4 + %0.2f CO2 --> %0.2f CH3OH + H2O\n\n',...
        solution.x(co2_idx)/solution.x(ch4_idx),...
        -solution.x(meoh_idx)/solution.x(ch4_idx),...
        -solution.x(h2o_idx)/solution.x(ch4_idx))

%Print the yield coefficient (grams biomass per mole CH4 produced)
fprintf('Predicted Yield Coefficient: %0.3f gDCW/mol CH4\n\n',-solution.f*1000/solution.x(ch4_idx))


% Print out the gibbs free energy prediction
% Add a warning for simulations with no concentrations given
if warning_flag
    warning('All external metabolite concentrations set to 1 mM');
end
fprintf('Predicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)

end
end