function model = createLatestModel()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file is meant to keep all model changes sequentially, starting at the
%beginning of May 2014

%Input: original model (M_maripaludis_methanogen_model.mat)
%Output: most current version of the model, with all changes chronicled


%First, load the model
load('M_maripaludis_methanogen_model.mat')

%Model variable is "M_mar"

%5/02/2014
%Add the 25 reactions additionally marked for addition
%4 Reactions from "other modified reactions"
model = addReaction(M_mar,'Nitrogen_fixation',...
    '16 H2O_c0 + 16 ATP_c0 + 8 H_e0 + 8 Reducedferredoxin_c0 + N2_c0 <=> 6 H_c0 + 16 Phosphate_c0 + 16 ADP_c0 + 2 NH3_c0 + 8 Oxidizedferredoxin_c0 + H2_c0');
model = addReaction(model,'ADP_glucokinase',...
    'ADP_c0 + D-Glucose_c0 <=> D-glucose-6-phosphate_c0 + AMP_c0');
model = addReaction(model,'PFK',...
    'ADP_c0 + D-fructose-6-phosphate_c0 -> D-fructose-1_6-bisphosphate_c0 + AMP_c0');
model = addReaction(model,'EC_1.2.1',...
    'CO2_c0 + Acetyl-CoA_c0 + NADH_c0 <=> CoA_c0 + NAD_c0 + Pyruvate_c0');

%21 Reactions from "added reactions"
%Named reactions first
model = addReaction(model,'Asparagine_biosynthesis_II',...
    'ATP_c0 + L-Aspartate_c0 + NH3_c0 <=> PPi_c0 + AMP_c0 + L-Asparagine_c0');
model = addReaction(model,'Glutamate_biosynthesis_III',...
    '2-Oxoglutarate_c0 + H_e0 + NADPH_c0 + NH3_c0 -> L-Glutamate_c0 + H2O_c0 + NADP_c0');
model = addReaction(model,'CoM_5',...
    'L-Cysteine_c0 + sulfoacetaldehyde_c0 -> H2O_c0 + 2-(sulfomethyl)thiazolidine-4-carboxylate_c0');
model = addReaction(model,'CoM_6',...
    '2-(sulfomethyl)thiazolidine-4-carboxylate_c0 -> sulfoethylcysteine_c0');
model = addReaction(model,'CoM_7',...
    'H2O_c0 + sulfoethylcysteine_c0 -> NH3_c0 + CoM_c0 + Pyruvate_c0');
model = addReaction(model,'Aconitase',...
    'Citrate_c0  -> H2O_c0 + cis-Aconitate_c0');
model = addReaction(model,'Aconitase_II',...
    'H2O_c0 + cis-Aconitate_c0  -> D-threo-Isocitrate_c0');
model = addReaction(model,'Isocitrate_dehydrogenase',...
    'NADP_c0 + D-threo-Isocitrate_c0  -> CO2_c0 + NADPH_c0 + 2-oxoglutarate_c0');
model = addReaction(model,'Formate_hydrogenlyase',...
    'Formate_c0 + H_e0  -> CO2_c0 + H2_c0');
model = addReaction(model,'H(2)-dependent methylenetetrahydromethanopterin dehydrogenase',...
    'H2_c0 + Reduced_coenzyme_F420_c0 + 5_10-Methenyltetrahydromethanopterin_c0  <=> H_e0 + Coenzyme_F420_c0 + 5_10-Methylenetetrahydromethanopterin_c0');
model = addReaction(model,'Citrate_synthase',...
    'Citrate_c0 + CoA_c0 + H_e0 <=> Acetyl-CoA_c0 + H2O_c0 + Oxaloacetate_c0');
model = addReaction(model,'Sulfopyruvate decarboxylase',...
    '3-sulfopyruvate_c0 + H_e0 -> CO2_c0 + sulfoacetaldehyde_c0');
model = addReaction(model,'Sulfolactate dehydrogenase',...
    '(R)-sulfolactate_c0 + NAD_c0 -> NADH_c0 + H_e0 + 3-sulfopyruvate_c0');
model = addReaction(model,'2_Phosphosulfolactate_phosphohydrolase',...
    '2R-Phosphosulfolactate_c0 + H2O_c0 -> Phosphate_c0 + H_c0 + (R)-sulfolactate_c0');
model = addReaction(model,'Diaminopimelate_aminotransferase',...
    'tetrahydrodipicolinate_c0 + H_e0 + H2O_c0 + L-Glutamate_c0 -> 2-oxoglutarate_c0 + LL-2_6-Diaminopimelate_c0');
model = addReaction(model,'Citramalate_synthase',...
    'Acetyl-CoA_c0 + H2O_c0 + Pyruvate_c0 -> H_e0 + CoA_c0 + Citramalate_c0');
model = addReaction(model,'Isopropylmalate_isomerase_I',...
    'Citramalate_c0 -> H2O_c0 + Citraconate_c0');
model = addReaction(model,'Isopropylmalate_isomerase_II',...
    'H2O_c0 + Citraconate_c0 -> beta-methyl-D-malate_c0');
model = addReaction(model,'Methylmalate_dehydrogenase',...
    'beta-methyl-D-malate_c0 + NAD_c0 -> NADH_c0 + CO2_c0 + 2-Oxobutyrate_c0');
model = addReaction(model,'Acetohydroxybutanoate_synthase',...
    '2-Oxobutyrate_c0 + H_e0 + Pyruvate_c0 -> 2-Aceto-2-hydroxybutanoate_c0 + CO2_c0');
model = addReaction(model,'Acetohydroxy_acid_isomeroreductase',...
    'NADPH_c0 + H_e0 + 2-Aceto-2-hydroxybutanoate_c0 -> 2_3-Dihydroxy-3-methylvalerate_c0 + NADP_c0');

%Associate genes with added reactions
model = changeGeneAssociation(model,'Asparagine_biosynthesis_II','mmp0918');
model = changeGeneAssociation(model,'Diaminopimelate_aminotransferase','mmp1527');
%Flag on this one: said mmp008...is that supposed to be mmp0081?
model = changeGeneAssociation(model,'Glutamate_biosynthesis_III','mmp0080 or mmp0082 or mmp0496');
model = changeGeneAssociation(model,'Citramalate_synthase','mmp1018');
model = changeGeneAssociation(model,'Isopropylmalate_isomerase_I','mmp1149');
model = changeGeneAssociation(model,'Isopropylmalate_isomerase_II','mmp1480');
model = changeGeneAssociation(model,'Methylmalate_dehydrogenase','mmp0539');
model = changeGeneAssociation(model,'Acetohydroxybutanoate_synthase','mmp0650 or mmp0651');
model = changeGeneAssociation(model,'Acetohydroxy_acid_isomeroreductase','mmp0654');
model = changeGeneAssociation(model,'2_Phosphosulfolactate_phosphohydrolase','mmp0161');
model = changeGeneAssociation(model,'Sulfolactate dehydrogenase','mmp1133');
model = changeGeneAssociation(model,'Sulfopyruvate decarboxylase','mmp0411 or mmp1689');
model = changeGeneAssociation(model,'Aconitase','mmp1480');
model = changeGeneAssociation(model,'Aconitase_II','mmp1480');
model = changeGeneAssociation(model,'Isocitrate_dehydrogenase','mmp0880');
model = changeGeneAssociation(model,'Formate_hydrogenlyase','mmp1298');
model = changeGeneAssociation(model,'H(2)-dependent methylenetetrahydromethanopterin dehydrogenase','mmp0127');
model = changeGeneAssociation(model,'Nitrogen_fixation','mmp0853 or mmp0856 or mmp0857');
model = changeGeneAssociation(model,'ADP_glucokinase','mmp1296');
model = changeGeneAssociation(model,'PFK','mmp1296');


%Remove 3 reactions from "other modified reactions"
model = removeRxns(model,{'rxn00216_c0','rxn03079_c0','rxn06874_c0'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End creation of initial model (05/02/2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
