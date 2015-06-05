function methaneFVA(model,num)

%Run FVA of 100%, see variability of methanogenesis

%First run FVA of 100%
[fva_min,fva_max] = fluxVariability(model,num);

%Pull out the indices of important things:
[~,h2_idx]  = intersect(model.rxns,'Ex_cpd11640_c0');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011_e0');
[~,ch4_idx] = intersect(model.rxns,'Ex_cpd01024_c0');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001_e0');
[~,form_idx] = intersect(model.rxns,'EX_cpd00047_e0');
[~,bio_idx] = intersect(model.rxnNames,'EX Biomass c0');

%Print out the flux range for each one 
fprintf('\n\nBiomass flux range: %0.2f-%0.2f\n',fva_min(bio_idx),fva_max(bio_idx));
fprintf('CO2 flux: %0.2f-%0.2f\n',fva_min(co2_idx),fva_max(co2_idx))
fprintf('H2 flux: %0.2f-%0.2f\n',fva_min(h2_idx),fva_max(h2_idx))
fprintf('H2O flux: %0.2f-%0.2f\n',fva_min(h2o_idx),fva_max(h2o_idx))
fprintf('CH4 flux: %0.2f-%0.2f\n',fva_min(ch4_idx),fva_max(ch4_idx))
fprintf('Formate flux range: %0.2f-%0.2f\n',fva_min(form_idx),fva_max(form_idx))