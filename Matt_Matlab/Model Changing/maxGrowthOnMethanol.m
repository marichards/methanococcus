function [solution,gibbs_flux,model] = maxGrowthOnMethanol(model,substrate_rxns,concentrations)

%Add in the reaction for converting methanol to methyl-coM
model = addReaction(model,{'rxn10568_c0','Methanol: coenzyme M methyltransferase'},...
    'Methanol_c0 + CoM_c0 <=> Methyl_CoM_c0 + H2O_c0');

%Add in the methanol uptake
model = addReaction(model,{'rxn10570_c0','Methanol diffusion'},...
    'Methanol_e0 <=> Methanol_c0');
model = addReaction(model,{'EX_cpd00116_e0','EX_Methanol_e0'},...
    'Methanol_e0 <=> ');

% Add metabolite info for methanol
[~,idx] = intersect(model.mets,'Methanol_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='CH4O';
[~,idx] = intersect(model.mets,'Methanol_e0');
model.metCharge(idx)=0;
model.metFormulas{idx}='CH4O';

% Add it to the free energy
[~,idx] = intersect(model.rxns,'EX_cpd00116_e0');
model.freeEnergy(idx) = -0.0302;

%Turn off the CO2 uptake
%model = changeRxnBounds(model,'EX_cpd00011_e0',0,'l');
% Turn on the acetate uptake
%model = changeRxnBounds(model,'EX_cpd00029_e0',-1000,'l');

% Turn off MFR too
model = changeRxnBounds(model,'rxn11938_c0',0,'b');

% Specify substrate reactions and concentrations as 1 mM if not given
if nargin<2
    
    substrate_rxns = {'EX_cpd00116_e0','EX_cpd11640_e0','EX_cpd01024_e0'};
    concentrations = [1 1 1];
end
%Solve by maximizing biomass
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001_e0');

%Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640_e0');
[~,meoh_idx] = intersect(model.rxns,'EX_cpd00116_e0');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011_e0');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024_e0');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001_e0');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013_e0');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009_e0');


%Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);

if solution.f > 0 
%Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))

%Print the per-CO2 actual reaction
fprintf('\nOverall reaction:\nCH3OH + H2 --> H2O + CH4\n')
fprintf('\nModel overall reaction (per mole CH4)\n')
fprintf('%0.2f CH3OH + %0.2f H2 --> %0.2f H2O + CH4\n\n',-solution.x(meoh_idx)/solution.x(ch4_idx),...
    -solution.x(h2_idx)/solution.x(ch4_idx),solution.x(h2o_idx)/solution.x(ch4_idx))

%Print the yield coefficient (grams biomass per mole CH4 produced)
fprintf('Predicted Yield Coefficient: %0.3f gDCW/mol CH4\n\n',solution.f*1000/solution.x(ch4_idx))

%Find the ATP reaction index
[~,atp_idx] = intersect(model.rxns,'ATPS');
%Print the ATP yield coefficient (ATP per CH4)
fprintf('Expected ATP Yield: 0.5\n')
fprintf('Predicted ATP Yield: %0.3f\n\n', solution.x(atp_idx)/solution.x(ch4_idx))

%Print out the gibbs free energy prediction
fprintf('Predicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)

end
end