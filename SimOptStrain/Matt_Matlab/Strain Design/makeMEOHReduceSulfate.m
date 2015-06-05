function solution = maxGrowthOnMethane(model)

%Turn off the MFR
%model = changeRxnBounds(model,'rxn11938_c0',0,'b');

%Overall methanol reaction(s): 
%H2O + SO4 + 5 CH4 = H2S + 5 CH3OH
%CH4 + 4 NO3- -> CO2 + 4 NO2- + 2 H2O

% From paper:
%CH4 + SO4(-2) -> HCO3- + HS- + H2O
%CH4 + 4 NO3- -> CO2 + 4 NO2- + 2 H2O

%Change bounds such that methane goes IN instead of OUT
model = changeRxnBounds(model,'Ex_cpd01024_c0',-1000,'l');
model = changeRxnBounds(model,'Ex_cpd01024_c0',0,'u');
%Change methanol to come OUT instead of IN
model = changeRxnBounds(model,'EX_methanol_c0',0,'l');
model = changeRxnBounds(model,'EX_methanol_c0',1000,'u');

%Turn off Hydrogen input and let it come out
model = changeRxnBounds(model,'Ex_cpd11640_c0',0,'l');
model = changeRxnBounds(model,'Ex_cpd11640_c0',1000,'u');

%Simulate growth
solution = optimizeCbModel(model,[],'one');

%Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'Ex_cpd11640_c0');
[~,meoh_idx] = intersect(model.rxns,'EX_methanol_c0');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011_e0');
[~,ch4_idx] = intersect(model.rxns,'Ex_cpd01024_c0');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001_e0');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013_e0');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009_e0');
[~,so4_idx] = intersect(model.rxns,'EX_cpd00048_e0');


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
fprintf('SO4 flux: %f\n',solution.x(so4_idx))

%Print the per-CO2 actual reaction
fprintf('\nExpected Overall reaction:\nH2O + SO4 + 5 CH4 = H2S + 5 CH3OH\n')
fprintf('\nModel overall reaction (per mole CH4)\n')
fprintf('%0.2f H2O + %0.2f SO4 + CH4 --> %0.2f H2S + CH3OH\n\n',-solution.x(h2o_idx)/solution.x(ch4_idx),...
    -solution.x(so4_idx)/solution.x(ch4_idx),...
    solution.x(h2s_idx)/solution.x(ch4_idx),...
    solution.x(meoh_idx)/solution.x(ch4_idx))

%Print the yield coefficient (grams biomass per mole CH4 produced)
fprintf('Predicted Yield Coefficient: %0.3f gDCW/mol CH4\n\n',solution.f*1000/solution.x(ch4_idx))

%Find the ATP reaction index
[~,atp_idx] = intersect(model.rxns,'ATPS');
%Print the ATP yield coefficient (ATP per CH4)
fprintf('Expected ATP Yield: 0.5\n')
fprintf('Predicted ATP Yield: %0.3f\n\n', solution.x(atp_idx)/solution.x(ch4_idx))

end