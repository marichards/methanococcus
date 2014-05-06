function solution = maxGrowthOnFormate(model)

%Simulate growth on CO2 and Formate media, print out the growth rate and
%relevant fluxes, return the full solution

%Note: Formate is HCO2
%Turn off H2 Input
model = changeRxnBounds(model,'Ex_cpd11640_c0',0,'l');
%Turn on Formate Input
model = changeRxnBounds(model,'EX_cpd00047_e0',-1000,'l');

%Find indices of important reactions
[~,h2_idx]  = intersect(model.rxns,'Ex_cpd11640_c0');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011_e0');
[~,ch4_idx] = intersect(model.rxns,'Ex_cpd01024_c0');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001_e0');
[~,form_idx] = intersect(model.rxns,'EX_cpd00047_e0');

%Solve for maximum biomass
solution = optimizeCbModel(model,[],'one');

%Print out the fluxes
%Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);
%Print the reaction fluxes
fprintf('Formate flux: %f\n',solution.x(form_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))

%Try max CH4
model = changeObjective(model,'Ex_cpd01024_c0');

%Repeat above

%Solve for maximum biomass
solution = optimizeCbModel(model,[],'one');

%Find biomass index
[~,bio_idx] = intersect(model.rxns,'EX Biomass c0');
%Print out the fluxes
fprintf('\nMaximize Methane Production')
%Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',bio_idx);
%Print the reaction fluxes
fprintf('Formate flux: %f\n',solution.x(form_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))