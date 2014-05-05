function solution = maxGrowthOnCO2(model)

%Simulate growth on CO2 and H2 media, print out the growth rate and
%relevant fluxes, return the full solution

%Solve by maximizing biomass
solution = optimizeCbModel(model,[],'one');

%Pull out the overall reaction CO2 + 4H2 --> CH4 + 2H2O

%Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'Ex_cpd11640_c0');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011_e0');
[~,ch4_idx] = intersect(model.rxns,'Ex_cpd01024_c0');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001_e0');


%Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);
%Print the reaction fluxes
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))

%Print the per-CO2 actual reaction
fprintf('\nOverall reaction:\nCO2 + 4 H2 --> 2 H2O + CH4\n')
fprintf('\nModel overall reaction (per mole CO2)\n')
fprintf('CO2 + %0.2f H2 --> %0.2f H2O + %0.2f CH4\n\n',solution.x(h2_idx)/solution.x(co2_idx),...
    -solution.x(h2o_idx)/solution.x(co2_idx),-solution.x(ch4_idx)/solution.x(co2_idx))

end