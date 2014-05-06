function solution = maxGrowthOnAcetate(model)

%Simulate growth on acetate media, print out the growth rate and
%relevant fluxes, return the full solution

%Turn off H2 Input
model = changeRxnBounds(model,'Ex_cpd11640_c0',0,'l');

%%%THIS CODE IS UNFINISHED BECAUSE THERE IS NO ACETATE UPTAKE RIGHT NOW

%Pull out the acetate reaction pathway...
%Acetate --> Acetyl-CoA (Acetate CoA Ligase)
[~,ligase_idx]=intersect(model.rxns,'rxn00175_c0');
%Acetyl-CoA --> Pyruvate
[~,pyr_idx]=intersect(model.rxns,'rxn05938_c0');
%Pyruvate --> Oxaloacetate
[~,oxac_idx]=intersect(model.rxns,'
