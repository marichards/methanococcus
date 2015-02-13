function solution = maxGrowthOnCO2(model)

%Simulate growth on CO2 and H2 media, print out the growth rate and
%relevant fluxes, return the full solution

%Turn down H2 to -45, not -1000
model = changeRxnBounds(model,'Ex_cpd11640_c0',-45,'l');
%Turn down CO2 to -12, not -1000
%model = changeRxnBounds(model,'EX_cpd00011_e0',-12,'l');

%Solve by maximizing biomass
solution = optimizeCbModel(model,[],'one');

%Pull out the overall reaction CO2 + 4H2 --> CH4 + 2H2O

%Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'Ex_cpd11640_c0');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011_e0');
[~,ch4_idx] = intersect(model.rxns,'Ex_cpd01024_c0');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001_e0');
[~,form_idx] = intersect(model.rxns,'EX_cpd00047_e0');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013_e0');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009_e0');
[~,ac_idx] = intersect(model.rxns,'EX_cpd00029_e0');
%[~,h2s_idx] = intersect(model.rxns,'EX_cpd00239_e0');


%Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);
%Print the reaction fluxes
fprintf('Formate flux: %f\n',solution.x(form_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))
fprintf('Acetate flux: %f\n',solution.x(ac_idx))
%fprintf('Sulfide flux: %f\n',solution.x(h2s_idx))


%Print the per-CO2 actual reaction
fprintf('\nOverall reaction:\nCO2 + 4 H2 --> 2 H2O + CH4\n')
fprintf('\nModel overall reaction (per mole CH4)\n')
fprintf('%0.2f CO2 + %0.2f H2 --> %0.2f H2O + CH4\n\n',-solution.x(co2_idx)/solution.x(ch4_idx),...
    -solution.x(h2_idx)/solution.x(ch4_idx),solution.x(h2o_idx)/solution.x(ch4_idx))

%Print the yield coefficient (grams biomass per mole CH4 produced)
fprintf('Measured Yield Coefficient: 2.31 +/- 0.29 gDCW/mol CH4\n')
fprintf('Predicted Yield Coefficient: %0.3f gDCW/mol CH4\n\n',solution.f*1000/solution.x(ch4_idx))

%Find the ATP reaction index
[~,atp_idx] = intersect(model.rxns,'ATPS');
%Print the ATP yield coefficient (ATP per CH4)
fprintf('Expected ATP Yield: 0.5\n')
fprintf('Predicted ATP Yield: %0.3f\n\n', solution.x(atp_idx)/solution.x(ch4_idx))

end