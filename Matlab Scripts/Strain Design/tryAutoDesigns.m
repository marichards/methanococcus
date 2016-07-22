function [solution_1,model_1,solution_2,model_2] = tryAutoDesigns(model)

% First, use below function to create methanol-producing path
model = addReverseMethanolPathway(model);

% Try the acetate transporter solution
[solution_1, ~, model_1] = addAcetateTransport(model);

% Try the phosphoglycerate solution
[solution_2, ~, model_2] = addPGPathway(model);

end


function [solution,gibbs_flux,model] = addAcetateTransport(model)

% Add the acetate transport and its metadata
model = addReaction(model,{'rxn05488[c0]',...
    'Acetate transport via proton symport'},...
    'cpd00029[c0] + cpd00067[c0] <=> cpd00029[e0] + cpd00067[e0]');
[~,idx] = intersect(model.rxns,'rxn05488[c0]');
model.freeEnergy(idx) = 0;
model = changeRxnBounds(model,'rxn05488[c0]',-inf,'l');
model = changeRxnBounds(model,'rxn05488[c0]',inf,'u');

% Specify substrate reactions and concentrations as 1 mM
substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]',...
    };
concentrations = [1 1 1];

% Solve by maximizing biomass, with negative free energy not allowed
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',false);

% Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(model.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');
[~,h2s_idx] = intersect(model.rxns,'EX_cpd00239[e0]');
[~,h_idx] = intersect(model.rxns,'EX_cpd00067[e0]');

if solution.f > 0 
% Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);

% Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))
fprintf('H2S flux: %f\n',solution.x(h2s_idx))
fprintf('H+ flux: %f\n',solution.x(h_idx))
end

% Print the Gibbs flux
fprintf('\nPredicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)

end

function [solution,gibbs_flux,model] = addPGPathway(model)

% Add the acetate transport and its metadata
model = addReaction(model,{'rxn00779[c0]',...
    'D-glyceraldehyde 3-phosphate:NADP+ oxidoreductase'},...
    'cpd00005[c0] + 2 cpd00067[c0] + cpd00169[c0] <=> cpd00001[c0] + cpd00006[c0] + cpd00102[c0]');
[~,idx] = intersect(model.rxns,'rxn00779[c0]');
model.freeEnergy(idx) = 0;
model = changeRxnBounds(model,'rxn00779[c0]',-inf,'l');
model = changeRxnBounds(model,'rxn00779[c0]',inf,'u');

% Specify substrate reactions and concentrations as 1 mM
substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]',...
    };
concentrations = [1 1 1];

% Solve by maximizing biomass, with negative free energy not allowed
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',false);

% Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(model.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');
[~,h2s_idx] = intersect(model.rxns,'EX_cpd00239[e0]');
[~,h_idx] = intersect(model.rxns,'EX_cpd00067[e0]');

if solution.f > 0 
% Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);

% Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))
fprintf('H2S flux: %f\n',solution.x(h2s_idx))
fprintf('H+ flux: %f\n',solution.x(h_idx))
end

% Print the Gibbs flux
fprintf('\nPredicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)

end

function model = addReverseMethanolPathway(model)

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
model = changeRxnBounds(model,'EX_cpd01024[e0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd01024[e0]',0,'u');
% Change methanol to come OUT instead of IN
% model = changeRxnBounds(model,'EX_cpd00116[e0]',10,'b');
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

% Let Eha/Ehb run at full blast
model = removeEhaBounds(model);

% Turn off Fwd, which will ruin everything if left on by allowing CO2 to be
% metabolized
model = changeRxnBounds(model,'rxn11938[c0]',0,'b');

end