function model = createLatestModel()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file is meant to keep all model changes sequentially, starting at the
%beginning of May 2014

%Input: original model (original_model.mat)
%Output: most current version of the model, with all changes chronicled


%First, load the model
load('original_model.mat')

%Model variable is "M_mar". Change it to "model"
model=M_mar;

%%%%%%%%%%%%%%%%%%
%5/09/2014
%%%%%%%%%%%%%%%%%%
%Remove the 'Unknown' and 'fig' genes
model=removeGene(model,'Unknown');
model=removeGene(model,'fig');

%Align all gene names with the convention of 'mmp#####'
%Do it in the genes
for i = 1:length(model.genes)
    %For each digit case, replace the "267377.1.peg." with "mmp" and the
    %correct number of 0s to make it 4 digits
    if regexp(model.genes{i},'267377.1.peg.[0-9]{1}$')
        model.genes{i} = regexprep(model.genes{i},'267377.1.peg.','mmp000');
    elseif regexp(model.genes{i},'267377.1.peg.[0-9]{2}$')
        model.genes{i} = regexprep(model.genes{i},'267377.1.peg.','mmp00');
    elseif regexp(model.genes{i},'267377.1.peg.[0-9]{3}$')
        model.genes{i} = regexprep(model.genes{i},'267377.1.peg.','mmp0');
    elseif regexp(model.genes{i},'267377.1.peg.[0-9]{4}$')
        model.genes{i} = regexprep(model.genes{i},'267377.1.peg.','mmp');
    end
end

%Do it in the rules
for i = 1:length(model.grRules)
    %For each digit case, replace the "267377.1.peg." with "mmp" and the
    %correct number of 0s to make it 4 digits
    if regexp(model.grRules{i},'\(*267377.1.peg.[0-9]{1}\)*$')
        model.grRules{i} = regexprep(model.grRules{i},'267377.1.peg.','mmp000');
    elseif regexp(model.grRules{i},'\(*267377.1.peg.[0-9]{2}\)*$')
        model.grRules{i} = regexprep(model.grRules{i},'267377.1.peg.','mmp00');
    elseif regexp(model.grRules{i},'\(*267377.1.peg.[0-9]{3}\)*$')
        model.grRules{i} = regexprep(model.grRules{i},'267377.1.peg.','mmp0');
    elseif regexp(model.grRules{i},'\(*267377.1.peg.[0-9]{4}\)*$')
        model.grRules{i} = regexprep(model.grRules{i},'267377.1.peg.','mmp');
    end
end

%%%%%%%%%%%%%%%%%%
%5/02/2014
%%%%%%%%%%%%%%%%%%
%Add the 25 reactions additionally marked for addition
%4 Reactions from "other modified reactions"
model = addReaction(model,'Nitrogen_fixation',...
    '16 H2O_c0 + 16 ATP_c0 + 8 H_e0 + 8 Reducedferredoxin_c0 + N2_c0 <=> 6 H_c0 + 16 Phosphate_c0 + 16 ADP_c0 + 2 NH3_c0 + 8 Oxidizedferredoxin_c0 + H2_c0');
model = addReaction(model,'ADP_glucokinase',...
    'ADP_c0 + D-Glucose_c0 <=> D-glucose-6-phosphate_c0 + AMP_c0');
model = addReaction(model,'PFK',...
    'ADP_c0 + D-fructose-6-phosphate_c0 -> D-fructose-1_6-bisphosphate_c0 + AMP_c0');
model = addReaction(model,'EC_1.2.1',...
    'CO2_c0 + Acetyl-CoA_c0 + NADH_c0 <=> CoA_c0 + NAD_c0 + Pyruvate_c0');

%21 Reactions from "added reactions"
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
model = changeGeneAssociation(model,'Glutamate_biosynthesis_III','mmp0080 or mmp0081 or mmp0082 or mmp0496');
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
%%%%%%%%%%%%%%%%%%
%05/13/2014
%%%%%%%%%%%%%%%%%%
%Change reaction directions (from metacyc(?))
%Rxn 11650
model = changeRxnBounds(model,'rxn11650_c0',-1000,'l');
model = changeRxnBounds(model,'rxn11650_c0',0,'u');
%Rxn 03535
model = changeRxnBounds(model,'rxn03535_c0',-1000,'l');
model = changeRxnBounds(model,'rxn03535_c0',0,'u');
%Rxn 03540
model = changeRxnBounds(model,'rxn03540_c0',-1000,'l');
model = changeRxnBounds(model,'rxn03540_c0',0,'u');

%%%%%%%%%%%%%%%%%%
%05/15/2014
%%%%%%%%%%%%%%%%%%
%Changes from Juan:
%5/13
%Add the reaction 'pyruvate-dependent arginine decarboxylase'
model = addReaction(model,'pyruvate-dependent arginine decarboxylase',...
    'L-Arginine_c0 + H_c0 -> CO2_c0 + Agmatine_c0 ');

%Associate the gene 'mmp1582' with reaction 'pyruvate-dependent arginine decarboxylase'
model = changeGeneAssociation(model,...
    'pyruvate-dependent arginine decarboxylase','mmp1582');

%5/14
%Associate the gene 'mmp1038 or mmp1039' with reaction 'ATP_synthase'
model = changeGeneAssociation(model,'ATP_synthase','mmp1038 or mmp1039');
%Associate the gene 'mmp0123' with reaction 'rxn03004_c0'
model = changeGeneAssociation(model,'rxn03004_c0','mmp0123');

%remove the gene association of 'mmp0882' with reaction 'rxn03084_c0'
%(gene mmp0178 or mmp0179 are still associated with this reaciton)
model = changeGeneAssociation(model,'rxn03004_c0','mmp0178 or mmp0179');

%Associate the gene 'mmp0882' with reaction 'rxn02937_c0'
%(gene mmp1254 are still associated with this reaction)
model = changeGeneAssociation(model,'rxn02937_c0','mmp1254 or mmp0882');

%5/15
%Add the reaction 'aldehyde dehydrogenase'
model = addReaction(model,'aldehyde dehydrogenase',...
    'Acetaldehyde_c0 + NAD_c0 + CoA_c0 -> H_c0 + NADH_c0 + Acetyl-CoA_c0 ');

%Associate the gene 'mmp1423' with reaction 'aldehyde dehydrogenase'
model = changeGeneAssociation(model,'aldehyde dehydrogenase','mmp1423');

%Associate the gene 'mmp0391 or mmp1527' with reaction 'rxn00260_c0'
model = changeGeneAssociation(model,'rxn00260_c0',...
'(mmp1396 or mmp1216 or mmp1072 or mmp0391 or mmp1527');

%Associate the gene '0082 or mmp0081 or mmp0080' with reaction 'rxn00085_c0'
model = changeGeneAssociation(model,'rxn00085_c0',...
    '(mmp0496 or mmp0082 or mmp0081 or mmp0080');

%Add gene mmp1259 and acossiated with rxn02269_c0
model = changeGeneAssociation(model,'rxn02269_c0','mmp1259');

%%%
%5/15 Changes from Me:
%Add CO-dehydrogenase
%Note: CO goes nowhere now, it is an orphan metabolite
%Should have a CO2+Fd_red-->CO+Fd_ox
model = addReaction(model,'CO dehydrogenase',...
    'CO_c0 + CoA_c0 + 5-Methyl-H4MPT_c0 -> H4MPT_c0 + Acetyl-CoA_c0');
%Associate it with mmp0980,0981,0983,0984,0985
model = changeGeneAssociation(model,'CO dehydrogenase',...
    'mmp0980 and mmp0981 and mmp0983 and mmp0984 and mmp0985');
%Add acetate exchange and transport, both reversible
model = addReaction(model,'EX_Acetate_e0',...
    'Acetate_e0 <=> ');
model = addReaction(model,'Acetate transport',...
    'Acetate_e0 <=> Acetate_c0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of 5/15 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



