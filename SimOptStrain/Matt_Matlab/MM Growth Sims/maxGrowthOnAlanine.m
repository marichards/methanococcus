function solution = maxGrowthOnAlanine(model)

%Simulate growth on CO2 and H2 media, print out the growth rate and
%relevant fluxes, return the full solution

%Turn down H2 to -10, not -1000
model = changeRxnBounds(model,'Ex_cpd11640_c0',-10,'l');
%Turn off NH3 input and turn on L-alanine (D is the same)
model = changeRxnBounds(model,'EX_cpd00013_e0',0,'l');
model = changeRxnBounds(model,'EX_cpd00035_e0',-10,'l');

%Solve by maximizing biomass
solution = optimizeCbModel(model,[],'one');

%Pull out the overall reaction CO2 + 4H2 --> CH4 + 2H2O

%Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'Ex_cpd11640_c0');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011_e0');
[~,ch4_idx] = intersect(model.rxns,'Ex_cpd01024_c0');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001_e0');
[~,form_idx] = intersect(model.rxns,'EX_cpd00047_e0');
[~,ala_idx] = intersect(model.rxns,'EX_cpd00035_e0');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009_e0');
[~,ac_idx] = intersect(model.rxns,'EX_cpd00029_e0');


%Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);
%Print the reaction fluxes
fprintf('Formate flux: %f\n',solution.x(form_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('Alanine flux: %f\n',solution.x(ala_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))
fprintf('Acetate flux: %f\n',solution.x(ac_idx))

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