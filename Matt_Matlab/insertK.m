function constrained = insertK(model,yield)

%Take in the M. maripaludis model, constrain it to match the desired yield,
%and return the constrained model. Note that yield is in gDCW/mol CH4

%Step 1: Find the biomass and methane reactions
[~,bio_idx] = intersect(model.rxns,'EX_cpd11416_c0');

%Change the metabolites of each reaction to add a "k" to the reactions
%Let's call it "Fake_met"
%Step 2: Modify our methane
constrained = addReaction(model,'Ex_cpd01024_c0','EX_Methane_c0 + Fake_met ->	');

%Step 3: Add it to the biomass using S matrix
[~,met_idx] = intersect(constrained.mets,'Fake_met');
constrained.S(met_idx,bio_idx) = 1000/yield;
