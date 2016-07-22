function [solution,gibbs_flux,model] =  produceMethanolReduceSulfite(model,substrate_rxns,concentrations)



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

% Change bounds to function in reverse
% Change bounds such that methane goes IN instead of OUT
model = changeRxnBounds(model,'EX_cpd01024[e0]',-10,'l');
model = changeRxnBounds(model,'EX_cpd01024[e0]',0,'u');
% Change methanol to come OUT instead of IN
model = changeRxnBounds(model,'EX_cpd00116[e0]',0,'l');
model = changeRxnBounds(model,'EX_cpd00116[e0]',1000,'u');
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

% Add Sulfite reduction to sulfide (F420 as carrier)
model = addReaction(model,'F420-dependent_sulfite_reductase',...
    'cpd00067[c0] + cpd00081[c0] + 3 cpd00792[c0]  -> 3 cpd00001[c0] + 3 cpd00649[c0] + cpd00239[c0]');

% Add Sulfite source, not worrying about exchange/transporter for now
model = addReaction(model,'Sulfite_supply','cpd00081[c0] <=> ');

% Give it pyruvate and let it GROW!
model = addReaction(model,'Pyruvate_supply','cpd00020[c0] <=> ');

% Add a row to freeEnergy to make dimensions correct
[~,idx] = intersect(model.rxns,'F420-dependent_sulfite_reductase');
model.freeEnergy(idx) = 0;
[~,idx] = intersect(model.rxns,'Sulfite_supply');
model.freeEnergy(idx) = 0;
[~,idx] = intersect(model.rxns,'Pyruvate_supply');
model.freeEnergy(idx) = 0;

% Specify substrate reactions and concentrations as 1 mM if not given
if nargin<2    
    substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]'};
    concentrations = [1 1 1];
    warning_flag = 1;
else
    warning_flag = 0;
end

% Solve by maximizing biomass
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]');

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

%Find the ATP reaction index
[~,atp_idx] = intersect(model.rxns,'ATPS');
%Print the ATP yield coefficient (ATP per CH4)
fprintf('Expected ATP Yield: 0.5\n')
fprintf('Predicted ATP Yield: %0.3f\n\n', solution.x(atp_idx)/solution.x(meoh_idx))

% Print out the gibbs free energy prediction
% Add a warning for simulations with no concentrations given
if warning_flag
    warning('All external metabolite concentrations set to 1 mM');
end
fprintf('Predicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)

end
end
