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
    model.genes{i} = regexprep(model.genes{i},'267377.1.peg.([0-9]{4})','mmp$1');
    model.genes{i} = regexprep(model.genes{i},'267377.1.peg.([0-9]{3})','mmp0$1');
    model.genes{i} = regexprep(model.genes{i},'267377.1.peg.([0-9]{2})','mmp00$1');
    model.genes{i} = regexprep(model.genes{i},'267377.1.peg.([0-9]{1})','mmp000$1');
    
    
    
end

%Do it in the rules
for i = 1:length(model.grRules)
    %For each digit case, replace the "267377.1.peg." with "mmp" and the
    %correct number of 0s to make it 4 digits
    model.grRules{i} = regexprep(model.grRules{i},'267377.1.peg.([0-9]{4})','mmp$1');
    model.grRules{i} = regexprep(model.grRules{i},'267377.1.peg.([0-9]{3})','mmp0$1');
    model.grRules{i} = regexprep(model.grRules{i},'267377.1.peg.([0-9]{2})','mmp00$1');
    model.grRules{i} = regexprep(model.grRules{i},'267377.1.peg.([0-9]{1})','mmp000$1');  
end

%%%%%%%%%%%%%%%%%%
%5/02/2014 (All changes but changing rxnNames)
%5/22/2014 (Change IDs to abbreviations, add rxnNames)
%%%%%%%%%%%%%%%%%%
%Add the 25 reactions additionally marked for addition
%2 Reactions from "other modified reactions"
%Modify the charge balance for the first one
model = addReaction(model,{'rxn06874_c0','reduced ferredoxin dinitrogen oxidoreductase ATP-hydrolysing c0'},...
    '16 H2O_c0 + 16 ATP_c0 + 4 Reducedferredoxin_c0 + N2_c0 <=> 16 Phosphate_c0 + 16 ADP_c0 + 2 NH3_c0 + 4 Oxidizedferredoxin_c0 + H2_c0 + 6 H_c0');
model = addReaction(model,{'rxn04042_c0','ADP:D-glucose 6-phosphotransferase'},...
    'ADP_c0 + D-Glucose_c0 <=> D-glucose-6-phosphate_c0 + AMP_c0');
model = addReaction(model,{'rxn04043_c0','ADP:D-fructose-6-phosphate 1-phosphotransferase'},...
    'ADP_c0 + D-fructose-6-phosphate_c0 <=> D-fructose-1_6-bisphosphate_c0 + AMP_c0');

%21 Reactions from "added reactions"
model = addReaction(model,{'rxn00340_c0','L-Aspartate:ammonia ligase (AMP-forming)'},...
    'ATP_c0 + L-Aspartate_c0 + NH3_c0 -> PPi_c0 + AMP_c0 + L-Asparagine_c0');
model = addReaction(model,{'rxn00184_c0','L-Glutamate:NADP+ oxidoreductase (deaminating)'},...
    '2-Oxoglutarate_c0 + H_c0 + NADPH_c0 + NH3_c0 <=> L-Glutamate_c0 + H2O_c0 + NADP_c0');
model = addReaction(model,{'rxn10603_c0','Sulfomethyl thiazolidine synthase'},...
    'L-Cysteine_c0 + sulfoacetaldehyde_c0 <=> H2O_c0 + 2-(sulfomethyl)thiazolidine-4-carboxylate_c0');
model = addReaction(model,{'rxn10598_c0','Sulfoethylcysteine synthase'},...
    '2-(sulfomethyl)thiazolidine-4-carboxylate_c0 + NADH_c0 + H_c0 <=> sulfoethylcysteine_c0 + NAD_c0');
model = addReaction(model,{'rxn10479_c0','Coenzyme M synthase'},...
    'H2O_c0 + sulfoethylcysteine_c0 <=> NH3_c0 + CoM_c0 + Pyruvate_c0');
model = addReaction(model,{'rxn00974_c0','Citrate hydrolase'},...
    'Citrate_c0  <=> H2O_c0 + cis-Aconitate_c0');
model = addReaction(model,{'rxn01388_c0','Isocitrate hydrolase'},...
    'H2O_c0 + cis-Aconitate_c0  <=> D-threo-Isocitrate_c0');
model = addReaction(model,{'rxn00198_c0','Isocitrate:NADP+ oxidoreductase (decarboxylating)'},...
    'NADP_c0 + D-threo-Isocitrate_c0  <=> CO2_c0 + NADPH_c0 + 2-Oxoglutarate_c0');
model = addReaction(model,{'rxn08518_c0','Formate-hydrogen lyase'},...
    'Formate_c0 + H_c0  -> CO2_c0 + H2_c0');
model = addReaction(model,{'rxn03079_c0','F420:5,10-methenyltetrahydromethanopterin oxidoreductase'},...
    'H_c0 + Coenzyme_F420_c0 + 5_10-Methylenetetrahydromethanopterin_c0  <=> Reduced_coenzyme_F420_c0 + 5_10-Methenyltetrahydromethanopterin_c0');
model = addReaction(model,{'rxn00256_c0','Citrate oxaloacetate-lyase'},...
    'Citrate_c0 + CoA_c0 + H_c0 <=> Acetyl-CoA_c0 + H2O_c0 + Oxaloacetate_c0');
model = addReaction(model,{'rxn04026_c0','Sulfopyruvate carboxy-lyase'},...
    '3-sulfopyruvate_c0 + H_c0 -> CO2_c0 + sulfoacetaldehyde_c0');
%Modified (R)-sulfolactate to (2R)-3-sulfolactate to match the seed on
%5/20/2014
model = addReaction(model,{'rxn04934_c0','(R)-2-hydroxyacid:NAD+ oxidoreductase'},...
    '(2R)-3-sulfolactate_c0 + NAD_c0 <=> NADH_c0 + H_c0 + 3-sulfopyruvate_c0');
%Modified (R)-sulfolactate to (2R)-3-sulfolactate to match the seed on
%5/20/2014
model = addReaction(model,{'rxn04036_c0','(R)-2-phospho-3-sulfolactate phosphohydrolase'},...
    '2R-Phosphosulfolactate_c0 + H2O_c0 -> Phosphate_c0 + H_c0 + (2R)-3-sulfolactate_c0');
model = addReaction(model,{'rxn07741_c0','LL-2,6-diaminoheptanedioate:2-Oxoglutarate aminotransferase'},...
    'tetrahydrodipicolinate_c0 + H_c0 + H2O_c0 + L-Glutamate_c0 <=> 2-Oxoglutarate_c0 + LL-2_6-Diaminopimelate_c0');
model = addReaction(model,{'rxn05109_c0','Citramalate_synthase'},...
    'Acetyl-CoA_c0 + H2O_c0 + Pyruvate_c0 -> H_c0 + CoA_c0 + Citramalate_c0');
model = addReaction(model,{'rxn02749_c0','(R)-2-Methylmalate hydro-lyase'},...
    'Citramalate_c0 <=> H2O_c0 + Citraconate_c0');
%In following 2 reactions, replaced beta-methyl-d-malate with D-erythro-3-methylmalate
model = addReaction(model,{'rxn02751_c0','Isopropylmalate_isomerase'},...
    'H2O_c0 + Citraconate_c0 + 2 H_c0 <=> D-erythro-3-methylmalate_c0 + H2_c0');
model = addReaction(model,{'rxn00735_c0','Methylmalate_dehydrogenase'},...
    'D-erythro-3-methylmalate_c0 + NAD_c0 + H2_c0 <=> NADH_c0 + CO2_c0 + 2-Oxobutyrate_c0 + 2 H_c0');
model = addReaction(model,{'rxn08043_c0','2-aceto-2-hydroxybutanoate synthase'},...
    '2-Oxobutyrate_c0 + H_c0 + Pyruvate_c0 -> 2-Aceto-2-hydroxybutanoate_c0 + CO2_c0');
model = addReaction(model,{'rxn08764_c0','ketol-acid reductoisomerase (2-Acetolactate)'},...
    'NADPH_c0 + H_c0 + 2-Aceto-2-hydroxybutanoate_c0 <=> 2_3-Dihydroxy-3-methylvalerate_c0 + NADP_c0');

%Associate genes with added reactions
model = changeGeneAssociation(model,'rxn06874_c0','mmp0853 and mmp0856 and mmp0857 and mmp0858 and mmp0859 and mmp0860');
model = changeGeneAssociation(model,'rxn04042_c0','mmp1296');
model = changeGeneAssociation(model,'rxn04043_c0','mmp1296');

model = changeGeneAssociation(model,'rxn00340_c0','mmp0918');
model = changeGeneAssociation(model,'rxn00184_c0','mmp0080 or mmp0081 or mmp0082 or mmp0496');



model = changeGeneAssociation(model,'rxn00974_c0','mmp1480');
model = changeGeneAssociation(model,'rxn01388_c0','mmp1480');
model = changeGeneAssociation(model,'rxn00198_c0','mmp0880');
model = changeGeneAssociation(model,'rxn08518_c0','mmp1298');
model = changeGeneAssociation(model,'rxn03079_c0','mmp0127');

model = changeGeneAssociation(model,'rxn04026_c0','mmp0411 or mmp1689');
model = changeGeneAssociation(model,'rxn04934_c0','mmp1133');
model = changeGeneAssociation(model,'rxn04036_c0','mmp0161');
model = changeGeneAssociation(model,'rxn07741_c0','mmp1527');
model = changeGeneAssociation(model,'rxn05109_c0','mmp1018');
model = changeGeneAssociation(model,'rxn02749_c0','mmp1149');
model = changeGeneAssociation(model,'rxn02751_c0','mmp1480');
model = changeGeneAssociation(model,'rxn00735_c0','mmp0539');
model = changeGeneAssociation(model,'rxn08043_c0','mmp0650 or mmp0651');
model = changeGeneAssociation(model,'rxn08764_c0','mmp0654');

%Remove reaction from "other modified reactions"
model = removeRxns(model,{'rxn00216_c0'});

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
model = addReaction(model,{'rxn00405_c0','Arginine carboxy-lyase'},...
    'L-Arginine_c0 + H_c0 -> CO2_c0 + Agmatine_c0 ');
%Associate the gene 'mmp1582' with reaction 'pyruvate-dependent arginine decarboxylase'
model = changeGeneAssociation(model,...
    'rxn00405_c0','mmp1582');

%5/14
%Associate the gene 'mmp1038 or mmp1039' with reaction 'ATP_synthase'
%Note: rxnID will change later, but right here it's still 'ATP_synthase'
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
model = addReaction(model,{'rxn00171_c0','Acetaldehyde:NAD+ oxidoreductase (CoA-acetylating'},...
    'Acetaldehyde_c0 + NAD_c0 + CoA_c0 <=> H_c0 + NADH_c0 + Acetyl-CoA_c0 ');
%Associate the gene 'mmp1423' with reaction 'aldehyde dehydrogenase'
model = changeGeneAssociation(model,'rxn00171_c0','mmp1423');

%Associate the gene 'mmp0391 or mmp1527' with reaction 'rxn00260_c0'
model = changeGeneAssociation(model,'rxn00260_c0',...
'(mmp1396 or mmp1216 or mmp1072 or mmp0391 or mmp1527)');

%Associate the gene '0082 or mmp0081 or mmp0080' with reaction 'rxn00085_c0'
model = changeGeneAssociation(model,'rxn00085_c0',...
    '(mmp0496 or mmp0082 or mmp0081 or mmp0080)');

%Add gene mmp1259 and acossiated with rxn02269_c0
model = changeGeneAssociation(model,'rxn02269_c0','mmp1259');

%%%
%5/15 Changes from Me:
%Add CO-dehydrogenase
model = addReaction(model,{'CODH_ACS','Carbon monoxide dehydrogenase'},...
    'CO2_c0 + CoA_c0 + 2 H_c0 + Reducedferredoxin_c0 + 5-Methyl-H4MPT_c0 <=> Acetyl-CoA_c0 + Oxidizedferredoxin_c0 + H2O_c0 + H4MPT_c0');
%Associate it with mmpmmp0980,0981,0983,0984,0985
model = changeGeneAssociation(model,'CODH_ACS',...
    'mmp0979 and mmp0980 and mmp0981 and mmp0982 and mmp0983 and mmp0984 and mmp0985');
%Add acetate exchange and transport, both reversible

%%9/17 change: turn exchange off for now!
model = addReaction(model,{'EX_cpd00029_e0','Acetate exchange'},...
    'Acetate_e0 <=> ');
model = changeRxnBounds(model,'EX_cpd00029_e0',0,'l');
model = addReaction(model,{'ACT','Acetate transport'},...
    'Acetate_e0 <=> Acetate_c0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of 5/15 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5/20 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add charges for added metabolites
%There are 13 added metabolites...add charges 1 by 1 (from Kbase)
%    'sulfoacetaldehyde_c0'
[~,idx] = intersect(model.mets,'sulfoacetaldehyde_c0');
model.metCharge(idx)=-1;
%    '2-(sulfomethyl)thiazolidine-4-carboxylate_c0'
[~,idx] = intersect(model.mets,'2-(sulfomethyl)thiazolidine-4-carboxylate_c0');
model.metCharge(idx)=-1;
%    'sulfoethylcysteine_c0'
[~,idx] = intersect(model.mets,'sulfoethylcysteine_c0');
model.metCharge(idx)=-1;
%    'Citrate_c0'
[~,idx] = intersect(model.mets,'Citrate_c0');
model.metCharge(idx)=-3;
%    'cis-Aconitate_c0'
[~,idx] = intersect(model.mets,'cis-Aconitate_c0');
model.metCharge(idx)=-3;
%    'D-threo-Isocitrate_c0'
[~,idx] = intersect(model.mets,'D-threo-Isocitrate_c0');
model.metCharge(idx)=-3;
%    '3-sulfopyruvate_c0'
[~,idx] = intersect(model.mets,'3-sulfopyruvate_c0');
model.metCharge(idx)=-2;
%    '(2R)-3-sulfolactate_c0'
[~,idx] = intersect(model.mets,'(2R)-3-sulfolactate_c0');
model.metCharge(idx)=-2;
%    'Citramalate_c0'
[~,idx] = intersect(model.mets,'Citramalate_c0');
model.metCharge(idx)=-2;
%    'Citraconate_c0'
[~,idx] = intersect(model.mets,'Citraconate_c0');
model.metCharge(idx)=-2;
%    'D-erythro-3-Methylmalate_c0'
[~,idx] = intersect(model.mets,'D-erythro-3-Methylmalate_c0');
model.metCharge(idx)=-2;
%    'Acetate_e0'
[~,idx] = intersect(model.mets,'Acetate_e0');
model.metCharge(idx)=-1;

%Fix charges for 3 other reactions
%Reaction 07191_c0; change 2 Fd to 1
%Find each rxn name first
[~,idx] = intersect(model.rxns,'rxn07191_c0');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn07191_c0',name},...
    'H2O_c0 + Glyceraldehyde3-phosphate_c0 + Oxidizedferredoxin_c0 <=> 3.000000 H_c0 + 3-Phosphoglycerate_c0 + Reducedferredoxin_c0');
%Reaction 04045_c0; Remove 2 protons from the left and add H2
[~,idx] = intersect(model.rxns,'rxn04045_c0');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn04045_c0',name},...
    'Sirohydrochlorin_c0 + Co2_c0 + H2_c0	<=>	Cobalt-precorrin_2_c0');
%Reaction 05029_c0; add 2 protons to the right and H2 to the left
[~,idx] = intersect(model.rxns,'rxn05029_c0');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn05029_c0',name},...
    'ATP_c0 + Cobinamide_c0 + H2_c0 ->	Triphosphate_c0 + Adenosyl_cobinamide_c0 + 2 H_c0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5/23 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find the ATP synthase reaction and change its ID
[~,idx] = intersect(model.rxnNames,'ATP_synthase');
model.rxns{idx} = 'rxn08173_c0';

%%%%
%Juan's changes from 5/16

%add gene mmp0341 and acossiated with rxn00250_c0
model = changeGeneAssociation(model,'rxn00250_c0','mmp0341 or mmp0340');

%add gene mmp0049 and acossiated with rxn00102_c0
model = changeGeneAssociation(model,'rxn00102_c0','mmp1299 or mmp0049');

%change the gene association of rxn00371_c0 from mmp1233 to (mmp0138 and mmp0139) or (mmp1297 and mmp1298)
model = changeGeneAssociation(model,'rxn00371_c0','(mmp0138 and mmp0139) or (mmp1297 and mmp1298)');

%add gene mmp1274 and acossiated with rxn00175_c0
model = changeGeneAssociation(model,'rxn00175_c0','mmp0148 or mmp1274');

%change all the 'H_e0' to 'H_c0' in all the added reactions. 
%In this case, it is interesting to find that f score changed, from 4.21 to 3.17
%,maybe indicating the effects to flux. 
%%%NOTE: Matt already did this! %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5/30 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Change the media
%Take out reactions
%On 7/24, don't remove, turn them off
%EX Urea e0	-1000
model = changeRxnBounds(model,'EX_cpd00073_e0',0,'l');
%EX fe3 e0	-1000
%model = changeRxnBounds(model,'EX_cpd10516_e0',0,'l');
%EX Spermine e0	-1000
model = changeRxnBounds(model,'EX_cpd00558_e0',0,'l');
%EX Nitrate e0	-1000
model = changeRxnBounds(model,'EX_cpd00209_e0',0,'l');
%EX BET e0	-1000
model = changeRxnBounds(model,'EX_cpd00540_e0',0,'l');
%EX Cytosine e0	-1000
model = changeRxnBounds(model,'EX_cpd00307_e0',0,'l');
%EX glycogenn-1 c0	-1000
model = changeRxnBounds(model,'EX_cpd15302_c0',0,'l');
%EX Uracil e0	-1000
model = changeRxnBounds(model,'EX_cpd00092_e0',0,'l');
%EX ddca e0	-1000
model = changeRxnBounds(model,'EX_cpd01741_e0',0,'l');
%EX Dephospho-CoA e0	-1000
model = changeRxnBounds(model,'EX_cpd00655_e0',0,'l');
%EX Cobinamide e0	-1000 %Necessary for Colamide?
model = changeRxnBounds(model,'EX_cpd03422_e0',0,'l');
%EX Oxidized glutathione e0	-1000
model = changeRxnBounds(model,'EX_cpd00111_e0',0,'l');

%Metals in the biomass that must be removed
%For now, don't remove their transporters or exchanges
%EX Zn2 e0	-1000
%model = changeRxnBounds(model,'EX_cpd00034_e0',0,'l');
%EX Cu2 e0	-1000
%model = changeRxnBounds(model,'EX_cpd00058_e0',0,'l');
%EX Mn2 e0	-1000
%model = changeRxnBounds(model,'EX_cpd00030_e0',0,'l');

%Change the biomass to take these out
%Find these in the model
[~,idx] = intersect(model.mets,{'Zn2_c0','Cu2_c0','Mn2_c0'});
%Remove them from the biomass
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(idx,bio_idx)=0;

%Turn rxn00175 irreversible forward
model=changeRxnBounds(model,'rxn00175_c0',0,'l');
%Turn rxn00549 reversible to allow growth
model=changeRxnBounds(model,'rxn00549_c0',-1000,'l');

%Note: I haven't removed exchanges or transporters for the 3 metals yet!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6/18 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Change rxn00781_c0 from NAD to NADP (John's suggestion)
[~,idx] = intersect(model.rxns,'rxn00781_c0');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn00781_c0',name},...
    'Phosphate_c0 + NADP_c0 + Glyceraldehyde3-phosphate_c0 <=> NADPH_c0 + 1_3-Bisphospho-D-glycerate_c0');

%Change Gene Associations to AND for ATP Synthase and rxn00250_c0, plus add
%some genes to ATP Synthase
model = changeGeneAssociation(model,'rxn08173_c0',...
    'mmp1038 and mmp1039 and mmp1040 and mmp1041 and mmp1042 and mmp1043 and mmp1044 and mmp1045 and mmp1046');
model = changeGeneAssociation(model,'rxn00250_c0','mmp0340 and mmp0341');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6/23 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Turn off  rxn07191_c0, make rxn05398_c0 and CODH_ACS irrev (John spreadsheet 6/19)
%model = changeRxnBounds(model,'rxn07191_c0',0,'b');
%model = changeRxnBounds(model,'rxn05938_c0',0,'l');
%model = changeRxnBounds(model,'CODH_ACS',0,'l');

%Add genes for Eha: mmp1448-1467
%%%%%9/30: Change name to Eha/Ehb, add more genes
%Original statement commented
%model = changeGeneAssociation(model,'Eha',...
%    'mmp1448 and mmp1449 and mmp1450 and mmp1451 and mmp1452 and mmp1453
%    and mmp1454 and mmp1455 and mmp1456 and mmp1457 and mmp1458 and
%    mmp1459 and mmp1460 and mmp1461 and mmp1462 and mmp1463 and mmp1464 and mmp1465 and mmp1466 and mmp1467')
[~,idx] = intersect(model.rxns,'Eha');
model.rxns{idx} = 'Eha/Ehb';
model.rxnNames{idx} = 'Eha/Ehb';
model = changeGeneAssociation(model,'Eha/Ehb',...
    '(mmp1448 and mmp1449 and mmp1450 and mmp1451 and mmp1452 and mmp1453 and mmp1454 and mmp1455 and mmp1456 and mmp1457 and mmp1458 and mmp1459 and mmp1460 and mmp1461 and mmp1462 and mmp1463 and mmp1464 and mmp1465 and mmp1466 and mmp1467) or (mmp1621 and mmp1622 and mmp1623 and mmp1624 and mmp1625 and mmp1626 and mmp1627 and mmp1628 and mmp1629 and mmp1073 and mmp1074 and mmp1469 and mmp0400)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6/26 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Group EhA, HdrABC, and rxn11938_c0 together by giving them a specific
%ferredoxin


%Add EhB, indolepyruvate oxidoreductase
%%%%%
%9/30: Don't add Ehb, just change genes and make Eha both of them
%Original statement
%model = addReaction(model,'Ehb',...
%    'Oxidizedferredoxin_c0 + 2.000000 Na_e0 + H2_c0 <=>	2.000000 H_c0 + Reducedferredoxin_c0 + 2.000000 Na_c0');
model = addReaction(model,{'rxn10561_c0','Indolepyruvate ferredoxin oxidoreductase'},...
    'Indole-3-pyruvate_c0 + Oxidizedferredoxin_c0 + CoA_c0 <=> S-2-(indol-3-yl)acetyl-CoA_c0 + CO2_c0 + Reducedferredoxin_c0 + H_c0');

%Add Genes for HdrABC and take them from rxn03126_c0
model = changeGeneAssociation(model,'HdrABC',...
    '(mmp0825 or mmp1697) and ((mmp1053 and mmp1054) or (mmp1155 and mmp1154))');
model = changeGeneAssociation(model,'rxn03126_c0','');

%Add Genes for EhB and indolepyruvate oxidoreductase
%model = changeGeneAssociation(model,'Ehb',...
%    'mmp1621 and mmp1622 and mmp1623 and mmp1624 and mmp1625 and mmp1626 and mmp1627 and mmp1628 and mmp1629 and mmp1073 and mmp1074 and mmp1469 and mmp0400');
model = changeGeneAssociation(model,'rxn10561_c0',...
    '(mmp0315 and mmp0316) or (mmp0713 and mmp0714)');

%%%%%%%%%%%%%%%%%%%
%8/5/2014 Changes
%%%%%%%%%%%%%%%%%%%
%Change mmp0372 to the rxn03079_c0
model = changeGeneAssociation(model,'rxn03079_c0','mmp0372');
%Make a new reaction for mmp0127 and associate it
model = addReaction(model,{'rxn06696_c0','Hydrogen:5,10-methenyltetrahydromethanopterin oxidoreductase'},...
    '5_10-Methenyltetrahydromethanopterin_c0 + H2_c0 <=> 5_10-Methylenetetrahydromethanopterin_c0 + H_c0');
model = changeGeneAssociation(model,'rxn06696_c0','mmp0127');
%Change name of rxn00371_c0 to Formate F420 oxidoreductase
[~,idx] = intersect(model.rxns,'rxn00371_c0');
model.rxnNames{idx} = 'Formate F420 oxidoreductase';

%%%%%%%%%%%%%%%%%%%
%8/8/2014 Changes
%%%%%%%%%%%%%%%%%%%
%Add E-matrix
model = addMetFormulas(model);

%%%%%%%%%%%%%%%%%%%
%8/19/2014 Changes
%%%%%%%%%%%%%%%%%%%
%Add charges for 2 from IPOR
%    'Indole-3-pyruvate_c0'
[~,idx] = intersect(model.mets,'Indole-3-pyruvate_c0');
model.metCharge(idx)=-1;
%    'S-2-(indol-3-yl)acetyl-CoA_c0'
[~,idx] = intersect(model.mets,'S-2-(indol-3-yl)acetyl-CoA_c0');
model.metCharge(idx)=-3;

%%%%%%%%%%%%%%%%%%%
%8/27/2014 Changes
%%%%%%%%%%%%%%%%%%%

%Take out methanophenazine reaction, it doesn't belong here
model = removeRxns(model,'rxn03126_c0');
    
    
%%%%%%%%%%%%%%%%%%%
%8/27/2014 Changes
%%%%%%%%%%%%%%%%%%%   

%Add in the 27 genes the way I THINK they should be, just confirm it later
%Add the Vhc and Vhu genes to Hdr
%model = changeGeneAssociation(model,'HdrABC',...
%    '((mmp0825 or mmp1697) and ((mmp1053 and mmp1054) or (mmp1155 and mmp1154))) and ((mmp0821 and mmp0822 and mmp0823 and mmp0824) or (mmp1692 and mmp1693 and mmp1694 and mmp1695 and mmp1696))');
%Add mmp1556-1557 to rxn03127_c0 
%model = changeGeneAssociation(model,'rxn03127_c0',...
%    'mmp1555 or mmp1556 or mmp1557 or mmp1558 or mmp1559');
%
%%%%%%%%%%%%%%%%%%%
%This section was hasty...don't do it now
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%9/4/2014 Changes
%%%%%%%%%%%%%%%%%%%   

%Make changes from John's spreadsheet corrections
%Catabolic
%Change genes of rxn11938_c0
model = changeGeneAssociation(model,'rxn11938_c0','(mmp1244 and mmp1245 and mmp1246 and mmp1247 and mmp1248 and (mmp1249 or mmp0070)) or (mmp0508 and mmp0509 and mmp0510 and (mmp0511 or mmp0512))');
%Change genes of rxn03020_c0
model = changeGeneAssociation(model,'rxn03020_c0','mmp1560 and mmp1561 and mmp1562 and mmp1563 and mmp1564 and mmp1565 and mmp1566 and mmp1567');
%Change genes of rxn03127_c0
model = changeGeneAssociation(model,'rxn03127_c0','(mmp1555 and mmp1559 and mmp1558)');
%Change genes of HdrABC
model = changeGeneAssociation(model,'HdrABC','((mmp0825 and mmp0821 and mmp0822 and mmp0823 and mmp0824) or (mmp1697 and mmp1692 and mmp1693 and mmp1694 and mmp1695 and mmp1696)) and ((mmp1053 and mmp1054) or (mmp1155 and mmp1154))');
%Change genes of rxn06299_c0
model = changeGeneAssociation(model,'rxn06299_c0','(mmp1382 and mmp1383 and mmp1384 and mmp1385) or (mmp0817 and mmp0818 and mmp0819 and mmp0820)');
%Changed Eha and Ehb to reversible earlier 
%%%%%%%%
%Change formula for rxn03079_c0 to take out H2...may imbalance the mass
[~,idx] = intersect(model.rxns,'rxn03079_c0');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn03079_c0',name},...
'H_c0 + Coenzyme_F420_c0 + 5_10-Methylenetetrahydromethanopterin_c0  <=> Reduced_coenzyme_F420_c0 + 5_10-Methenyltetrahydromethanopterin_c0');
%Change ATP synthase from H+ to Na+ by removing the reaction 8173 and
%putting in my own Na+ one
model = removeRxns(model,'rxn08173_c0');
model = addReaction(model,{'ATPS','ATP synthase:Sodium utilizing'},...
    'Phosphate_c0 + ADP_c0 + 4 Na_e0 + H_c0 <=> 4 Na_c0 + H2O_c0 + ATP_c0');
%Give it the gene associations from the previous synthase
model = changeGeneAssociation(model,'ATPS',...
    'mmp1038 and mmp1039 and mmp1040 and mmp1041 and mmp1042 and mmp1043 and mmp1044 and mmp1045 and mmp1046');
%Change Na+ and H+ antiporter stoich to 2 protons and change genes
[~,idx] = intersect(model.rxns,'rxn05209_c0');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn05209_c0',name},...
    '2 H_c0 + Na_e0  <=> 2 H_e0 + Na_c0');
model = changeGeneAssociation(model,'rxn05938_c0',...
    'mmp0864 or mmp0679 or mmp0587 or mmp0707');
%Remove the rxn00371_c0 (Formate w/NAD)
%model = removeRxns(model,'rxn00371_c0');
%Add the Formate dehydrogenase 
model = addReaction(model,{'Fdh','Formate dehydrogenase'},...
    'Formate_c0 + Coenzyme_F420_c0 + H_c0 <=> CO2_c0 + Reduced_coenzyme_F420_c0');
model = changeGeneAssociation(model,'Fdh','(mmp0138 and mmp0139) or (mmp1297 and mmp1298)');

%Add the Formate:Hdr
model = addReaction(model,{'Hdr_formate','Formate-utilizing heterodisulfide reductase'},...
    'CoM-S-S-CoB_c0 + 2 Formate_c0 + Oxidizedferredoxin_c0 -> 2 CO2_c0 + CoM_c0 + HTP_c0 + Reducedferredoxin_c0');
%%%%
%Altered 12/16/2014 to include mmp1696 or mmp0821, the VhuD/VhcD genes
%%%%
model = changeGeneAssociation(model,'Hdr_formate','((mmp0825 or mmp1697) and ((mmp1053 and mmp1054) or (mmp1155 and mmp1154))) and ((mmp0138 and mmp0139) or (mmp1297 and mmp1298)) and (mmp1696 or mmp0821)');

%Biosynthesis reactions
%Add Malate reaction (NADP version)
model = addReaction(model,{'rxn00249_c0','(S)-Malate:NADP+ oxidoreductase'},...
    'L-Malate_c0 + NADP_c0  <=> H_c0 + Oxaloacetate_c0 + NADPH_c0');
%Add Genes for it
model = changeGeneAssociation(model,'rxn00249_c0','mmp0645');
%Add genes for rxn0175 and rxn05938
model = changeGeneAssociation(model,'rxn00175_c0','mmp0148');
model = changeGeneAssociation(model,'rxn05938_c0',...
    '(mmp1504 and mmp1502 and mmp1505 and mmp1507 and mmp1503 and mmp1506)');
%Remove rxn 00288 and replace with the fumarate reaction
%model = removeRxns(model,'rxn00288_c0');
model = addReaction(model,{'Tfr','Thiol:fumarate reductase'},...
    'Fumarate_c0 + CoM_c0 + HTP_c0  -> Succinate_c0 + CoM-S-S-CoB_c0');
%Add genes for Tfr
model = changeGeneAssociation(model,'Tfr','mmp1277 and mmp1067');
%Change genes for 05939 and 00085
model = changeGeneAssociation(model,'rxn05939_c0','mmp1687 and mmp1316 and mmp0003 and mmp1315');
model = changeGeneAssociation(model,'rxn00085_c0','mmp0082 and mmp0081 and mmp0080');

%Making it work: turn off ability to bring in protons from outside that can
%artificially pump sodium and make ATP
%model = changeRxnBounds(model,'EX_cpd00067_e0',0,'l')
%Turning this off doesn't make things happy....it would ruin formate!

%%%%%%%%%%%%%%%
%9/11 Changes
%%%%%%%%%%%%%%%

%Add alanine exchanges but don't let it uptake alanine for now
model = addReaction(model,{'EX_cpd00035_e0','EX L-Alanine e0'},...
    'L-Alanine_e0 	->	');
model = addReaction(model,{'EX_cpd00117_e0','EX D-Alanine e0'},...
    'D-Alanine_e0 	->	');

%Add alanine permeases for both, gene mmp1511
model = addReaction(model,{'rxn05215_c0','L-Alanine permease'},...
    'L-Alanine_e0 + Na_e0 <=> L-Alanine_c0 + Na_c0');
model = addReaction(model,{'rxn13660_c0','D-Alanine permease'},...
    'D-Alanine_e0 + Na_e0 <=> D-Alanine_c0 + Na_c0');

model = changeGeneAssociation(model,'rxn05215_c0','mmp1511');
model = changeGeneAssociation(model,'rxn13660_c0','mmp1511');

%Remove the glutamate dehydrogenase and force flux through GOGAT cycle
model = removeRxns(model,'rxn00184_c0');


%%%%%%%%%%%%%
%9/30/2014
%%%%%%%%%%%%%
%Add nitrogen diffusion at 0 lb for now
model = addReaction(model,{'EX_cpd00528_e0','EX Nitrogen e0'},...
    'N2_e0 -> ');
%Important: also turn the nitrogen fixation irreversible, lest it blow up
%and let ATP be formed massively:
model = changeRxnBounds(model,'rxn06874_c0',0,'l');
%%%%%%%%%%%%%
%1/21/2015
%%%%%%%%%%%%%
%Add a diffusion transporter
model = addReaction(model,{'rxn10577_c0','Nitrogen exchange, diffusion'},...
    'N2_e0 <=> N2_c0');

%Match pairs of NAD/NADP reactions with directions John suggested
%Malate synthesis
model = changeRxnBounds(model,{'rxn00248_c0','rxn00249_c0'},0,'u');
%Prephenate pair is already correct
%Homoserine synthesis
model = changeRxnBounds(model,{'rxn01301_c0','rxn01302_c0'},0,'u');
%Take out the terin reactions
model = removeRxns(model,{'rxn01313_c0','rxn01314_c0'});
%Take out the aspartate reactions and replace them with aspartate dehydrogenase
%Don't put in the dehydrogenase
%(rxns 04952 and 04953)
%model = addReaction(model,{'rxn04952_c0','L-aspartate:NAD+ oxidoreductase (deaminating)'},...
%    'H2O_c0 + NAD_c0 + L-Aspartate_c0 <=> Oxaloacetate_c0 + NH3_c0 + NADH_c0 + H_c0');
%model = addReaction(model,{'rxn04953_c0','L-aspartate:NADP+ oxidoreductase (deaminating)'},...
%    'H2O_c0 + NADP_c0 + L-Aspartate_c0 <=> Oxaloacetate_c0 + NH3_c0 + NADPH_c0 + H_c0');
%model = changeRxnBounds(model,{'rxn04952_c0','rxn04953_c0'},0,'u');
%model = changeGeneAssociation(model,'rxn04952_c0','mmp0737');
%model = changeGeneAssociation(model,'rxn04953_c0','mmp0737');
model = removeRxns(model,{'rxn05117_c0','rxn05119_c0'});

%%%%%%%%%%%%
%11/19/2014
%%%%%%%%%%%%
%Change the name of rxn05938_c0 to "Pyruvate Oxidoreductase"
[~,idx] = intersect(model.rxns,'rxn05938_c0');
model.rxnNames{idx} = 'Pyruvate Oxidoreductase';

%%%%%%%%%%%%%%%
%11/20/2014
%%%%%%%%%%%%%%%
%Already changed CODH-->CODH_ACS
%Change comDE to an AND relationship
model = changeGeneAssociation(model,'rxn04026_c0','mmp0411 and mmp1689');

%%%%%%%%%%%%%%%
%11/24/2014
%%%%%%%%%%%%%%%
%Add in a biomass value for coenzyme M (0.0030965)
[~,coM_idx] = intersect(model.mets,'CoM_c0');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(coM_idx,bio_idx) = -0.0030965;

%%%%%%%%%%%%%%%
%12/16/2014
%%%%%%%%%%%%%%%
%Take out the incorrect FDH that uses NAD/NADH
model = removeRxns(model,'rxn00371_c0');



%%%%%%%%%%%%%
%1/13/2015
%%%%%%%%%%%%%
%Add a maintenance factor and change the H2 so both match M.barkeri model
model = changeRxnBounds(model,'Ex_cpd11640_c0',-45,'l');
model = changeRxnBounds(model,'rxn00062_c0',2,'b');

%%%%%%%%%%%%%
%1/29/2015
%%%%%%%%%%%%%
%Remove the sulfate uptake and transport
model = removeRxns(model,{'rxn05651_c0','EX_cpd00048_e0'});

%Replace it with H2S uptake and transport
model = addReaction(model,{'rxn10541_c0','H2S transport (diffusion)'},...
'H2S_e0 <=> H2S_c0');
model = addReaction(model,{'EX_cpd00239_e0','EX_H2S_e0'},...
'H2S_e0 <=>');

%Make model work by removing sulfate from the biomass
[~,so4_idx] = intersect(model.mets,'Sulfate_c0');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(so4_idx,bio_idx) = 0;

%Also, allow it to turn H2S into Sulfite
%model = changeRxnBounds(model,'rxn00623_c0',-1000,'l');

%Turn off the Hdr_Formate when not growing on formate
model = changeRxnBounds(model,'Hdr_formate',0,'b');

%Put in the reaction for H2S to Sulfite
model = addReaction(model,'Dsr-LP',...%'H2S_c0 -> Sulfite_c0');
    'H2S_c0 + 3 H2O_c0 + NAD_c0 -> Sulfite_c0 + NADH_c0 + 2 H2_c0 + 2 H_c0 ');
    
%%%%%%%%%%%%%
%02/10/2015
%%%%%%%%%%%%%

%Remove things from the biomass based on John's advice:
mets = {'Menaquinone_8_c0';...
    'Ubiquinone-8_c0';...
    'Bactoprenyl_diphosphate_c0';...
%    '2-Demethylmenaquinone_8_c0';...
    'Peptidoglycan_polymer_n_subunits_c0';...
    'Isoheptadecanoylcardiolipin_B._subtilis_c0';...
    'Stearoylcardiolipin_B._subtilis_c0';...
    'Diisoheptadecanoylphosphatidylethanolamine_c0';...
    'Diisoheptadecanoylphosphatidylglycerol_c0';...
    'Dianteisoheptadecanoylphosphatidylethanolamine_c0';...
    'Dianteisoheptadecanoylphosphatidylglycerol_c0';...
    '5-Methyltetrahydrofolate_c0';...
    'Tetrahydrofolate_c0';...
    'core_oligosaccharide_lipid_A_c0';...
    'Phosphatidylglycerol_dioctadecanoyl_c0';...
    'Anteisoheptadecanoylcardiolipin_B._subtilis_c0';...
    'phosphatidylethanolamine_dioctadecanoyl_c0';...
    '10-Formyltetrahydrofolate_c0';...
    'Peptidoglycan_polymer_n-1_subunits_c0';...
%Remove other electron carriers that lead up to these too
%   '2-Demethylmenaquinol_8_c0';...
   'Menaquinone_7_c0';...
   'Menaquinone_7_e0';...
   'Ubiquinol-8_c0';...
   'Menaquinol_8_c0';...
   'Heme_c0';...
    };
%rxns_to_remove = {};
for i=1:length(mets)
    %Remove it from the biomass
    [~,met_idx] = intersect(model.mets,mets{i});
    [~,bio_idx] = intersect(model.rxns,'biomass0');
    model.S(met_idx,bio_idx) = 0;
    
    %Search for reactions that include it
    rxns = findRxnsFromMets(model,mets{i});
    %Remove them
    model = removeRxns(model,rxns);
end
%Also grab their genes
%temp_genes = findGenesFromRxns(model,rxns_to_remove);
%Make them linear
%genes={};
%for i=1:length(temp_genes)
%    genes = [genes;temp_genes{i}];
%end
%Make them unique
%genes = unique(genes);
%%Also remove folate from the model

%Add in the 17 reactions for Coenzyme B (HTP) synthesis
%Step 1: 2-oxoglutarate + Acetyl-CoA ->trans-homoaconitate (aksA)
model = addReaction(model,{'rxn10434_c0','alpha-ketoglutarate:CoA ligase'},...
    '2-Oxoglutarate_c0 + Acetyl-CoA_c0 <=> CoA_c0 + H_c0 + trans-homoaconitate_c0');
model = changeGeneAssociation(model,'rxn10434_c0','mmp0153');
%Step 2: trans-homoaconitate <=> S-homocitrate
model = addReaction(model,{'rxn10608_c0','trans-homoaconitate hydrolase'},...
    'H2O_c0 + trans-homoaconitate_c0 <=> S-homocitrate_c0');
%Step 3: S-homocitrate <=> cis-homoaconitate
model = addReaction(model,{'rxn10599_c0','(S)-homocitrate dehydratase'},...
    'S-homocitrate_c0 <=> H2O_c0 + cis-Homoaconitate_c0');
%Step 4: cis-homoaconitate <=> threo-isohomocitrate (aksD/E)
model = addReaction(model,{'rxn10472_c0','cis-homoaconitate hydrolase'},...
    'H2O_c0 + cis-Homoaconitate_c0 <=> threo-isohomocitrate_c0');
model = changeGeneAssociation(model,'rxn10472_c0','mmp1480 and mmp0381');
%Step 5: threo-isohomocitrate + NAD -> alpha-ketoadipate + CO2 + NADH
%(aksF)
model = addReaction(model,{'rxn10612_c0','threo-isohomocitrate dehydrogenase'},...
    'NAD_c0 + threo-isohomocitrate_c0 <=> 2-oxohexanedioic_acid_c0 + CO2_c0 + NADH_c0');
model = changeGeneAssociation(model,'rxn10612_c0','mmp0880');
%Step 6: alpha-ketoadipate + Acetyl-CoA -> R-homocitrate (aksA)
model = addReaction(model,{'rxn10433_c0','alpha-ketoadipate:CoA ligase'},...
    '2-oxohexanedioic_acid_c0 + Acetyl-CoA_c0 + H2O_c0 <=> (R)-(homo)2citrate_c0 + CoA_c0 + H_c0');
model = changeGeneAssociation(model,'rxn10433_c0','mmp0153');
%Step 7: R-(homo)2citrate <=> cis-(homo)2aconitate (aksD/E)
model = addReaction(model,{'rxn10595_c0','(R)-(homo)2citrate dehydratase'},...
    '(R)-(homo)2citrate_c0 <=> H2O_c0 + cis-(homo)2aconitate_c0');
model = changeGeneAssociation(model,'rxn10595_c0','mmp1480 and mmp0381');
%Step 8: cis-(homo)2aconitate <=> threo-iso(homo)2citrate (aksD/E)
model = addReaction(model,{'rxn10468_c0','cis-(homo)2aconitate hydrolase'},...
    'H2O_c0 + cis-(homo)2aconitate_c0 <=> (-)threo-iso(homo)2citrate_c0');
model = changeGeneAssociation(model,'rxn10468_c0','mmp1480 and mmp0381');
%Step 9: threo-iso(homo)2citrate + NAD -> alpha-ketopimelate + CO2 + NADH
%(aksF)
model = addReaction(model,{'rxn10610_c0','threo-iso(homo)2citrate dehydrogenase'},...
    '(-)threo-iso(homo)2citrate_c0 + NAD_c0 <=> 2-oxoheptanedioic_acid_c0 + CO2_c0 + NADH_c0');
model = changeGeneAssociation(model,'rxn10610_c0','mmp0880');
%Step 10: alpha-ketopimelate + Acetyl-CoA -> R-(homo)3citrate (aksA)
model = addReaction(model,{'rxn10435_c0','alpha-ketopimelate:CoA ligase'},...
    '2-oxoheptanedioic_acid_c0 + Acetyl-CoA_c0 + H2O_c0 <=> (R)-(homo)3citrate_c0 + CoA_c0 + H_c0');
model = changeGeneAssociation(model,'rxn10435_c0','mmp0153');
%Step 11: R-(homo)3citrate <=> cis-(homo)3aconitate (aksD/E)
model = addReaction(model,{'rxn10596_c0','(R)-(homo)3citrate dehydratase'},...
    '(R)-(homo)3citrate_c0 <=> H2O_c0 + cis-(homo)3aconitate_c0');
model = changeGeneAssociation(model,'rxn10596_c0','mmp1480 and mmp0381');
%Step 12: cis-(homo)3aconitate <=> threo-iso(homo)3citrate (aksD/E)
model = addReaction(model,{'rxn10469_c0','cis-(homo)3aconitate hydrolase'},...
    'H2O_c0 + cis-(homo)3aconitate_c0 <=> (-)threo-iso(homo)3citrate_c0');
model = changeGeneAssociation(model,'rxn10469_c0','mmp1480 and mmp0381');
%Step 13: threo-iso(homo)3citrate + NAD -> alpha-ketosuberate + CO2 + NADH
%(aksF)
model = addReaction(model,{'rxn10611_c0','threo-iso(homo)3citrate dehydrogenase'},...
    '(-)threo-iso(homo)3citrate_c0 + NAD_c0 <=> 2-Oxosuberate_c0 + CO2_c0 + NADH_c0');
model = changeGeneAssociation(model,'rxn10611_c0','mmp0880');
%Step 14: alpha-ketosuberate -> 7-oxoheptanoic acid + CO2
model = addReaction(model,{'rxn11855_c0','rxn11855'},...
    '2-Oxosuberate_c0 + H_c0 <=> 7-oxoheptanoic_acid_c0 + CO2_c0');
%Step 15: 7-oxoheptanoic acid + Sulfur + 2e- -> 7-mercaptoheptanoic acid
model = addReaction(model,{'rxn10424_c0','7-mercaptoheptanoate synthase'},...
    '7-oxoheptanoic_acid_c0 + H_c0 + L-Cysteine_c0 + NADH_c0 <=> 7-mercaptoheptanoic_acid_c0 + L-Serine_c0 + NAD_c0');
%Step 16: 7-mercaptoheptanoic acid + threonine + ATP -> 7-mercaptoheptanoylthreonine + ADP
model = addReaction(model,{'rxn10425_c0','7-mercaptoheptanoylthreonine synthase'},...
    '7-mercaptoheptanoic_acid_c0 + ATP_c0 + L-Threonine_c0 <=> 7-mercaptoheptanoylthreonine_c0 + ADP_c0 + H_c0 + Phosphate_c0');
%Step 17: 7-mercaptoheptanoylthreonine + ATP -> Coenzyme B + ADP
model = addReaction(model,{'rxn10475_c0','Coenzyme B synthase'},...
    '7-mercaptoheptanoylthreonine_c0 + ATP_c0 <=> ADP_c0 + HTP_c0');

%Add Coenzyme B to the biomass
[~,coB_idx] = intersect(model.mets,'HTP_c0');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(coB_idx,bio_idx) = -0.0030965;

%%%%%%%%%%%%%
%02/24/2015
%%%%%%%%%%%%%
%Remove the 4 reactions that complete the TCA cycle (we added all 4 above)
model = removeRxns(model,{'rxn00974_c0','rxn00198_c0','rxn01388_c0','rxn00256_c0'});

%Add H4MPT Synthesis Pathway
%Step 1(MptA): GTP + 2 H2O -> Formate + PPi + 7,8-dihydronepterin 2' :3'-cyclicphosphate
model = addReaction(model,{'MptA','MptA'},...
    'GTP_c0 + H2O_c0 -> Formate_c0 + PPi_c0 + 7,8-dihydronepterin_2_3-cyclicphosphate_c0 + H_c0');
model = changeGeneAssociation(model,'MptA','mmp0034');
%Step 2(MptB): 7,8-dihydronepterin 2' :3'-cyclicphosphate + H2O -> Dihydroneopterin phosphate + H+
model = addReaction(model,{'rxn10490_c0','7,8-dihydronepterin 2'':3''-cyclicphosphate hydrolase'},...
    '7,8-dihydronepterin_2_3-cyclicphosphate_c0 + H2O_c0 -> Dihydroneopterin_phosphate_c0');
model = changeGeneAssociation(model,'rxn10490_c0','mmp0230');
%Step 3:  Dihydroneopterin phosphate + H2O <=> Dihydroneopterin + H+ + Ppi
model = addReaction(model,{'rxn03168_c0','Dihydroneopterin monophosphate dephosphorylase'},...
    'Dihydroneopterin_phosphate_c0 + H2O_c0 <=> Dihydroneopterin_c0 + H_c0 + Phosphate_c0');
%Step 4(MptD):Dihydroneopterin -> 6-hydroxymethyl-7,8-dihydropterin + glycolaldehyde
model = addReaction(model,{'rxn02504_c0','Dihydroneopterin aldolase'},...
    'Dihydroneopterin_c0 <=> 6-hydroxymethyl_dihydropterin_c0 + Glycolaldehyde_c0');
model = changeGeneAssociation(model,'rxn02504_c0','mmp0243');
%Step 5(MptE):6-hydroxymethyl-7,8-dihydropterin + ATP ->
%6-hydroxymethyl-7,8-dihydropterin diphosphate + PPi
model = addReaction(model,{'rxn02503_c0','ATP:2-amino-4-hydroxy-6-hydroxymethyl-7,8-dihydropteridine'},...
    '6-hydroxymethyl_dihydropterin_c0 + ATP_c0 <=> 2-Amino-4-hydroxy-6-hydroxymethyl-7-8-dihydropteridinediphosphate_c0 + AMP_c0');
model = changeGeneAssociation(model,'rxn02503_c0','mmp0579');
%Step 6(MptG): 4-aminobenzoate + PRPP -> beta-RFA-P
model = addReaction(model,{'rxn10446_c0','beta-ribofuranosylaminobenzene 5''-phosphate synthase'},...
    'ABEE_c0 + PRPP_c0 <=> 4-(B-D-ribofuranosyl)aminobenzene_5-phosphate_c0 + CO2_c0 + PPi_c0');
model = changeGeneAssociation(model,'rxn10446_c0','mmp0279');
%Step 7(MptH): 6-HMDP + beta-RFA-P <=> Intermediate + PPi
model = addReaction(model,{'rxn10491_c0',''},...
    '2-Amino-4-hydroxy-6-hydroxymethyl-7-8-dihydropteridinediphosphate_c0 + 4-(B-D-ribofuranosyl)aminobenzene_5-phosphate_c0 <=> 7,8-dihydropterin-6-ylmethyl-4-(B-D-ribofuranosyl)_aminobenzene_5-phosphate_c0 + PPi_c0');
%Step 8: Intermediate + 2 ATP + S-Hydroxyglutarate + 2 H+ + 2 NADH + PRPP +
%2 S-adenosyl-L-methionine -> H4MPT + 2 ADP + 2 PPi + Pi + 2 NAD + 2
%S-adenosyl-L-homocysteine
model = addReaction(model,{'H4MPTs','Tetrahydromethanopterin synthase'},...
    '7,8-dihydropterin-6-ylmethyl-4-(B-D-ribofuranosyl)_aminobenzene_5-phosphate_c0 + 2 ATP_c0 + (S)-2-Hydroxyglutarate_c0 + 2 NADH_c0 + PRPP_c0 + 2 S-Adenosyl-L-methionine_c0 + H2O_c0 <=> H4MPT_c0 + 2 ADP_c0 + 2 PPi_c0 + Phosphate_c0 + 2 NAD_c0 + 2 S-Adenosyl-homocysteine_c0');

%H4MPT Synthesis requires ABEE and 2-Hydroxyglutarate...add reactions for
%each
%Gapfill for 2-hydroxyglutarate
model = addReaction(model,{'rxn10432_c0','(S)-alpha-hydroxyglutarate dehydrogenase'},...
    'NADH_c0 + 2-Oxoglutarate_c0 + H_c0 <=> NAD_c0 + (S)-2-Hydroxyglutarate_c0');
%Path for ABEE synthesis (according to Porat et al)
%Need path to DKFP:
model = addReaction(model,{'rxn10494_c0','6-deoxy-5-ketofructose 1-phosphate synthase'},...
    'NADH_c0 + H_c0 + 4,5-diketo-6-deoxyfructose_1-phosphate_c0 <=> NAD_c0 + 6-deoxy-5-ketofructose-1-phosphate_c0');
model = addReaction(model,{'rxn10422_c0','4-ketofructose_1,6-biphosphate synthase'},...
    '4-ketofructose_1,6-bisphosphate_c0 <=> Phosphate_c0 + 4,5-diketo-6-deoxyfructose_1-phosphate_c0');
model = addReaction(model,{'rxn10520_c0','fructose-bisphoshate oxidase'},...
    'NAD_c0 + D-fructose-1_6-bisphosphate_c0 <=> 3 H_c0 + NADH_c0 + 4-ketofructose_1,6-bisphosphate_c0');
%From DKFP, add path to ABEE using genes
%aroA' gene: http://www.uniprot.org/uniprot/Q57843
model = addReaction(model,{'ADTHs','ADTH synthase'},...
    '6-deoxy-5-ketofructose-1-phosphate_c0 + L-Aspartate4-semialdehyde_c0 -> 2-amino-2,3,7-trideoxy-D-lyxo-hept-6-ulosonate_c0 + Hydroxypyruvaldehyde_phosphate_c0'); 
%Add aroA' gene from Porat 2006 paper:
model = changeGeneAssociation(model,'ADTHs','mmp0686');
%aroB' gene from Porat 2006 paper
model = addReaction(model,{'ADTHOR','ADTH oxidoreductase'},...
'2-amino-2,3,7-trideoxy-D-lyxo-hept-6-ulosonate_c0 + H2O_c0 + NAD_c0 + H_c0 <=> 3-Dehydroquinate_c0 + NH3_c0 + NADH_c0');
model = changeGeneAssociation(model,'ADTHOR','mmp0006');

%Add 5 hypothetical reactions with hypothetical compounds, all from Porat
%2006
%1) DHQ + NH3 -> 4-amino-DHQ {Hypothetical from Porat 2006}
model = addReaction(model,{'3DHQAT','3-dehydroquinate aminotransferase'},...
'3-Dehydroquinate_c0 + NH3_c0  <=> 4-Amino-3-Dehydroquinate_c0 + H_c0 + H2O_c0');
%2) 4-amino-DHQ -> 4-amino-DHS + H2O {Hypothetical from Porat 2006}
model = addReaction(model,{'4ADSs','4-aminodehydroshikimate synthase'},...
'4-Amino-3-Dehydroquinate_c0  <=> 4-Amino-3-Dehydroshikimate_c0 + H2O_c0');
%3) 4-amino-DHS + NADH -> 4-Aminoshikimate + NAD {Hypothetical from Porat 2006}
model = addReaction(model,{'4ASDH','4-aminoshikimate dehydrogenase'},...
'4-Amino-3-Dehydroshikimate_c0 + NADH_c0 + H_c0 <=>  4-Aminoshikimate_c0 + NAD_c0');
%4) 4-Aminoshikimate  -> 4-amino-3-hydroxycyclohexa-1,5-diene-1-carboxylate + H2O {Hypothetical from Porat 2006}
model = addReaction(model,{'4ASDHT','4-aminoshikimate dehydratase'},...
'4-Aminoshikimate_c0 <=> 4-amino-3-hydroxycyclohexa-1,5-diene-1-carboxylate_c0 + H2O_c0');
%5) 4-amino-3-hydroxycyclohexa-1,5-diene-1-carboxylate -> ABEE + H2O {Hypothetical from Porat 2006}
model = addReaction(model,{'ABEEs','4-aminobenzoic acid synthase'},...
'4-amino-3-hydroxycyclohexa-1,5-diene-1-carboxylate_c0 <=>  ABEE_c0 + H2O_c0');

%For now, add an artificial outlet for the Hydroxypyruvaldehyde_phosphate
model = addReaction(model,'SINK','Hydroxypyruvaldehyde_phosphate_c0 <=> ');

%Finally, add H4MPT to the biomass
[~,H4MPT_idx] = intersect(model.mets,'H4MPT_c0');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(H4MPT_idx,bio_idx) = -0.0030965;

%We need to consume the Glycolaldehyde_c0....manually GapFill it using the
%rxnlikelihoods file
model = addReaction(model,{'rxn00979_c0','Glycolaldehyde:NAD+ oxidoreductase'},...
    'H2O_c0 + NAD_c0 + Glycolaldehyde_c0 <=> NADH_c0 + 2 H_c0 + Glycolate_c0');
model = addReaction(model,{'rxn00512_c0','Glycolate:NAD+ oxidoreductase'},...
    'Glycolate_c0 + NAD_c0 <=> Glyoxalate_c0 + H_c0 + NADH_c0');
model = addReaction(model,{'rxn00272_c0','L-Alanine:glyoxylate aminotransferase'},...
    'L-Alanine_c0 + Glyoxalate_c0 <=> Glycine_c0 + Pyruvate_c0');

%Add formulas and charges for things
%L-Alanine_e0
[~,idx] = intersect(model.mets,'L-Alanine_e0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C3H7NO2';
%D-Alanine_e0
[~,idx] = intersect(model.mets,'D-Alanine_e0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C3H7NO2';
%N2_e0
[~,idx] = intersect(model.mets,'N2_e0');
model.metCharge(idx)=0;
model.metFormulas{idx}='N2';
%H2S_e0
[~,idx] = intersect(model.mets,'H2S_e0');
model.metCharge(idx)=0;
model.metFormulas{idx}='H2S';
%trans-homoaconitate_c0
[~,idx] = intersect(model.mets,'trans-homoaconitate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C7H5O6';
%S-homocitrate_c0
[~,idx] = intersect(model.mets,'S-homocitrate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C7H7O7';
%cis-Homoaconitate_c0
[~,idx] = intersect(model.mets,'cis-Homoaconitate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C7H5O6';
%threo-isohomocitrate_c0
[~,idx] = intersect(model.mets,'threo-isohomocitrate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C7H7O7';
%2-oxohexanedioic_acid_c0
[~,idx] = intersect(model.mets,'2-oxohexanedioic_acid_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C6H6O5';
%(R)-(homo)2citrate_c0
[~,idx] = intersect(model.mets,'(R)-(homo)2citrate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C8H9O7';
%cis-(homo)2aconitate_c0
[~,idx] = intersect(model.mets,'cis-(homo)2aconitate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C8H7O6';
%(-)threo-iso(homo)2citrate_c0
[~,idx] = intersect(model.mets,'(-)threo-iso(homo)2citrate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C8H9O7';
%2-oxoheptanedioic_acid_c0
[~,idx] = intersect(model.mets,'2-oxoheptanedioic_acid_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C7H8O5';
%(R)-(homo)3citrate_c0
[~,idx] = intersect(model.mets,'(R)-(homo)3citrate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C9H11O7';
%cis-(homo)3aconitate_c0
[~,idx] = intersect(model.mets,'cis-(homo)3aconitate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C9H9O6';
%(-)threo-iso(homo)3citrate_c0
[~,idx] = intersect(model.mets,'(-)threo-iso(homo)3citrate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C9H11O7';
%2-Oxosuberate_c0
[~,idx] = intersect(model.mets,'2-Oxosuberate_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C8H10O5';
%7-oxoheptanoic_acid_c0
[~,idx] = intersect(model.mets,'7-oxoheptanoic_acid_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C7H11O3';
%7-mercaptoheptanoic_acid_c0
[~,idx] = intersect(model.mets,'7-mercaptoheptanoic_acid_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C7H13O2S';
%7-mercaptoheptanoylthreonine_c0
[~,idx] = intersect(model.mets,'7-mercaptoheptanoylthreonine_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C11H20NO4S';
%7,8-dihydronepterin_2_3-cyclicphosphate_c0
[~,idx] = intersect(model.mets,'7,8-dihydronepterin_2_3-cyclicphosphate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C9H11N5O6P';
%Dihydroneopterin_phosphate_c0
[~,idx] = intersect(model.mets,'Dihydroneopterin_phosphate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C9H13N5O7P';
%Dihydroneopterin_c0
[~,idx] = intersect(model.mets,'Dihydroneopterin_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C9H13N5O4';
%6-hydroxymethyl_dihydropterin_c0
[~,idx] = intersect(model.mets,'6-hydroxymethyl_dihydropterin_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C7H9N5O2';
%Glycolaldehyde_c0
[~,idx] = intersect(model.mets,'Glycolaldehyde_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C2H4O2';
%2-Amino-4-hydroxy-6-hydroxymethyl-7-8-dihydropteridinediphosphate_c0
[~,idx] = intersect(model.mets,'2-Amino-4-hydroxy-6-hydroxymethyl-7-8-dihydropteridinediphosphate_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C7H9N5O8P2';
%ABEE_c0
[~,idx] = intersect(model.mets,'ABEE_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C7H6NO2';
%4-(B-D-ribofuranosyl)aminobenzene_5-phosphate_c0
[~,idx] = intersect(model.mets,'4-(B-D-ribofuranosyl)aminobenzene_5-phosphate_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C11H14NO7P';
%7,8-dihydropterin-6-ylmethyl-4-(B-D-ribofuranosyl)_aminobenzene_5-phosphate_c0
[~,idx] = intersect(model.mets,'7,8-dihydropterin-6-ylmethyl-4-(B-D-ribofuranosyl)_aminobenzene_5-phosphate_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C18H21N6O8P';
%(S)-2-Hydroxyglutarate_c0
[~,idx] = intersect(model.mets,'(S)-2-Hydroxyglutarate_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C5H6O5';
%4,5-diketo-6-deoxyfructose_1-phosphate_c0
[~,idx] = intersect(model.mets,'4,5-diketo-6-deoxyfructose_1-phosphate_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C6H7O8P';
%6-deoxy-5-ketofructose-1-phosphate_c0
[~,idx] = intersect(model.mets,'6-deoxy-5-ketofructose-1-phosphate_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C6H9O8P';
%4-ketofructose_1,6-bisphosphate_c0
[~,idx] = intersect(model.mets,'4-ketofructose_1,6-bisphosphate_c0');
model.metCharge(idx)=-4;
model.metFormulas{idx}='C6H8O12P2';
%2-amino-2,3,7-trideoxy-D-lyxo-hept-6-ulosonate_c0
[~,idx] = intersect(model.mets,'2-amino-2,3,7-trideoxy-D-lyxo-hept-6-ulosonate_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C7H11NO5';
%Hydroxypyruvaldehyde_phosphate_c0
[~,idx] = intersect(model.mets,'Hydroxypyruvaldehyde_phosphate_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C3H5O6P';
%3-Dehydroquinate_c0
[~,idx] = intersect(model.mets,'3-Dehydroquinate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C7H9O6';
%4-Amino-3-Dehydroquinate_c0
[~,idx] = intersect(model.mets,'4-Amino-3-Dehydroquinate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C7H10O5N';
%4-Amino-3-Dehydroshikimate_c0
[~,idx] = intersect(model.mets,'4-Amino-3-Dehydroshikimate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C7H8O4N'; 
%4-Aminoshikimate_c0
[~,idx] = intersect(model.mets,'4-Aminoshikimate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C7H10O4N'; 
%4-amino-3-hydroxycyclohexa-1,5-diene-1-carboxylate_c0
[~,idx] = intersect(model.mets,'4-amino-3-hydroxycyclohexa-1,5-diene-1-carboxylate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C7H8O3N'; 
%Glycolate_c0
[~,idx] = intersect(model.mets,'Glycolate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C2H3O3';
%Glyoxalate_c0
[~,idx] = intersect(model.mets,'Glyoxalate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C2HO3';

%%%%%%%%%%%%%
%3/2/2015
%%%%%%%%%%%%%
%Correct the final step of CoM synthesis
%Remove the incorrect reactions (from MetaCyc)
model = removeRxns(model,{'rxn10479_c0','rxn10598_c0','rxn10603_c0'});
%Add the hypothetical reaction from Graham et al
model = addReaction(model,{'COMs','Coenzyme M synthesis'},...
    'sulfoacetaldehyde_c0 + H2S_c0 + H2_c0 -> CoM_c0 + H2O_c0');
    
%Add steps for Methanofuran synthesis
%mfnA reaction
model = addReaction(model,{'rxn00529_c0','L-Tyrosine carboxylase'},...
'L-Tyrosine_c0 + H_c0 -> Tyramine_c0 + CO2_c0');
model = changeGeneAssociation(model,'rxn00529_c0','mmp0131');
%mfnD reaction
model = addReaction(model,{'MfnD','Tyramine-glutamate ligase'},...
'Tyramine_c0 + L-Glutamate_c0 + ATP_c0 -> gamma-Glutamyl-tyramine_c0 + ADP_c0 + Phosphate_c0 + H_c0');
%mfnB reaction
model = addReaction(model,{'MfnB','Tyramine-glutamate ligase'},...
'Glyceraldehyde3-phosphate_c0 + Glycerone-phosphate_c0 -> 4-(hydroxymethyl)-2-furancarboxaldehyde-phosphate_c0 + 2 H2O_c0 + Phosphate_c0 + H_c0');
model = changeGeneAssociation(model,'MfnB','mmp0708');
%mfnC reaction
model = addReaction(model,{'MfnC','4-HCF-P:alanine aminotransferase'},...
    '4-(hydroxymethyl)-2-furancarboxaldehyde-phosphate_c0 + L-Alanine_c0 -> 5-(aminomethyl)-3-furanmethanol-phosphate_c0 + Pyruvate_c0');
%3 Hypotheteical Reactions
model = addReaction(model,{'F1Pp','5-(aminomethyl)-3-furanmethanol-phosphate phosphorylase'},...
    '5-(aminomethyl)-3-furanmethanol-phosphate_c0 + ATP_c0 -> 5-(aminomethyl)-3-furanmethanol-pyrophosphate_c0 + ADP_c0');
model = addReaction(model,{'F1PPc','5-(aminomethyl)-3-furanmethanol-pyrophosphate condensation'},...
    '5-(aminomethyl)-3-furanmethanol-pyrophosphate_c0 + gamma-Glutamyl-tyramine_c0 -> 4((4-(2-aminoethyl)phenoxy)methyl)-2-furanmethanamine_c0 + PPi_c0 ');
model = addReaction(model,{'MFs','Methanofuran synthase'},...
    '4((4-(2-aminoethyl)phenoxy)methyl)-2-furanmethanamine_c0 + 3 L-Glutamate_c0 + H2_c0 -> Methanofuran_c0 + 2 NH3_c0 + 2 H2O_c0');

%Add Methanofuran to biomass
[~,MF_idx] = intersect(model.mets,'Methanofuran_c0');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(MF_idx,bio_idx) = -0.0030965;

%Need to add charges and formulae for new things
[~,idx] = intersect(model.mets,'Tyramine_c0');
model.metCharge(idx)=1;
model.metFormulas{idx}='C8H12NO';
[~,idx] = intersect(model.mets,'gamma-Glutamyl-tyramine_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C13H18N2O4';
[~,idx] = intersect(model.mets,'4-(hydroxymethyl)-2-furancarboxaldehyde-phosphate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C6H6O6P';
[~,idx] = intersect(model.mets,'5-(aminomethyl)-3-furanmethanol-phosphate_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C6H10O5NP';
[~,idx] = intersect(model.mets,'5-(aminomethyl)-3-furanmethanol-pyrophosphate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C6H10O8NP2';
[~,idx] = intersect(model.mets,'4((4-(2-aminoethyl)phenoxy)methyl)-2-furanmethanamine_c0');
model.metCharge(idx)=1;
model.metFormulas{idx}='C19H26O5N3';

%%%%%%%%%%%%%
%3/5/2015
%%%%%%%%%%%%%
%Add reactions for F430 synthesis
model = addReaction(model,{'rxn10508_c0','coenzyme F430 precursor synthase 1 (aminase)'},...
    '2 ATP_c0 + 2 H2O_c0 + 2 L-Glutamine_c0 + Ni2_c0 + Precorrin_2_c0 <=> 2 ADP_c0 + 4 H_c0 + 2 L-Glutamate_c0 + 2 Phosphate_c0 + Pyrrochorphinate_c0');
model = addReaction(model,{'rxn10509_c0','coenzyme F430 precursor synthase 2 (reduction)'},...
    '2 H_c0 + NADH_c0 + Pyrrochorphinate_c0 <=> NAD_c0 + Dihydrocorphinate_c0');
model = addReaction(model,{'rxn10510_c0','coenzyme F430 precursor synthase 3 (reduction)'},...
    'H_c0 + NADH_c0 + Dihydrocorphinate_c0 <=> NAD_c0 + Tetrahydrocorphinate_c0');
model = addReaction(model,{'rxn10511_c0','coenzyme F430 precursor synthase 4 (cyclization)'},...
    'Tetrahydrocorphinate_c0 <=> 15_17-seco-F430-17-acid_c0');
model = addReaction(model,{'rxn10512_c0','coenzyme F430 precursor synthase 5 (cyclization)'},...
    '15_17-seco-F430-17-acid_c0 + H_c0 <=> Factor_430_c0 + H2O_c0');
%Add Nickel to the media
model = addReaction(model,{'EX_cpd00244_e0','EX_Ni2_e0'},...
    'Ni2_e0 <=> ');
model = addReaction(model,{'rxn05174_c0','Nickel-ABC Transport'},...
    'H2O_c0 + ATP_c0 + Ni2_e0 -> ADP_c0 + Phosphate_c0 + H_c0 + Ni2_c0');

[~,idx] = intersect(model.mets,'Pyrrochorphinate_c0');
model.metCharge(idx)=-6;
model.metFormulas{idx}='C42H42N6NiO14';
[~,idx] = intersect(model.mets,'Dihydrocorphinate_c0');
model.metCharge(idx)=-5;
model.metFormulas{idx}='C42H45N6NiO14';
[~,idx] = intersect(model.mets,'Tetrahydrocorphinate_c0');
model.metCharge(idx)=-5;
model.metFormulas{idx}='C42H47N6NiO14';
[~,idx] = intersect(model.mets,'15_17-seco-F430-17-acid_c0');
model.metCharge(idx)=-5;
model.metFormulas{idx}='C42H47N6NiO14';
[~,idx] = intersect(model.mets,'Factor_430_c0');
model.metCharge(idx)=-4;
model.metFormulas{idx}='C42H46N6O13Ni';
[~,idx] = intersect(model.mets,'Ni2_c0');
model.metCharge(idx)=2;
model.metFormulas{idx}='Ni';
[~,idx] = intersect(model.mets,'Ni2_e0');
model.metCharge(idx)=2;
model.metFormulas{idx}='Ni';
%Add F430 to the biomass
[~,f430_idx] = intersect(model.mets,'Factor_430_c0');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(f430_idx,bio_idx) = -0.0030965;

%%%%%%%%%%%%%
%3/23/2015
%%%%%%%%%%%%%
%Turn off the ability to import protons and sodium ions from media
model = changeRxnBounds(model,'EX_cpd00067_e0',0,'l');
model = changeRxnBounds(model,'EX_cpd00971_e0',0,'l');

%Turn off the ability to dump out formate
model = changeRxnBounds(model,'EX_cpd00047_e0',0,'b');

%%%%%%%%%%%%%
%3/24/2015
%%%%%%%%%%%%%
%Add adenylate kinase
model = addReaction(model,{'rxn00097_c0','ATP:AMP phosphotransferase'},...
    'ATP_c0 + AMP_c0 -> 2 ADP_c0');
model = changeGeneAssociation(model,'rxn00097_c0','mmp1031');

%%%%%%%%%%%%%
%3/25/2015
%%%%%%%%%%%%%
%Turn off all alanine outlets
model = changeRxnBounds(model,'EX_cpd00035_e0',0,'u');
model = changeRxnBounds(model,'EX_cpd00117_e0',0,'u');

%%%%%%%%%%%%%
%3/26/2015
%%%%%%%%%%%%%
%Add in the F420 synthesis pathway
%Start from FO synthase(CofG/H)
model = addReaction(model,{'rxn10499_c0','7,8-didemethyl-8-hydroxy-5-deazariboflavin synthetase'},...
    'H2O_c0 + 2 NADP_c0 + p-hydroxyphenylpyruvate_c0 + 4-1-D-Ribitylamino-5-aminouracil_c0 <=> 2 NADPH_c0 + NH3_c0 + 2 H_c0 + Oxalate_c0 + 7,8-didemethyl-8-hydroxy-5-deazariboflavin_c0');
model = changeGeneAssociation(model,'rxn10499_c0','mmp0876 and mmp0056');
%Next add a pathway to make lactate, the lactate dehydrogenase
model = addReaction(model,{'rxn01053_c0','(S)-lactaldehyde:NAD+ oxidoreductase'},...
    'H2O_c0 + NAD_c0 + L-Lactaldehyde_c0 <=> NADH_c0 + 2 H_c0 + L-Lactate_c0');
model = changeGeneAssociation(model,'rxn01053_c0','mmp1487');
%We already get L-Lactaldehyde_c0 from L-fuculose-1-phosphate in the model
%(mmp1187)
%Next add LLPG synthesis (CofC)
model = addReaction(model,{'rxn10567_c0','LPPG synthetase'},...
    'GTP_c0 + H_c0 + 2-phospho-L-lactate_c0 <=> PPi_c0 + lactyl-(2)-diphospho-(5)-guanosine_c0');
model = changeGeneAssociation(model,'rxn10567_c0','mmp0117');
%This requires 2-phospholactate; kinase is unknown, but it's there. Add a
%hypothetical reaction (a manual gapfill of sorts), the only one in Kbase:
model = addReaction(model,{'rxn10420_c0','2-phospho-L-lactate synthase'},...
    'GTP_c0 + L-Lactate_c0 <=> 2-phospho-L-lactate_c0 + GDP_c0 + H_c0');
%Add LLPG + FO reaction, the F420-0 synthase (CofD):
model = addReaction(model,{'rxn10566_c0','LPPG:Fo 2-phospho-L-lactate transferase'},...
    '7,8-didemethyl-8-hydroxy-5-deazariboflavin_c0 + lactyl-(2)-diphospho-(5)-guanosine_c0 <=> GMP_c0 + Coenzyme_F420-0_c0');
model = changeGeneAssociation(model,'rxn10566_c0','mmp0404');
%Add glutamates to the F420-0 to make F420-2, the functional type (CofE)
model = addReaction(model,{'rxn10525_c0','gamma-F420-0:gamma-L-glutamate ligase'},...
    'L-Glutamate_c0 + GTP_c0 + Coenzyme_F420-0_c0 <=> Phosphate_c0 + GDP_c0 + H_c0 + Coenzyme_F420-1_c0');
model = changeGeneAssociation(model,'rxn10525_c0','mmp0937');
model = addReaction(model,{'rxn10526_c0','gamma-F420-1:gamma-L-glutamate ligase'},...
    'L-Glutamate_c0 + GTP_c0 + Coenzyme_F420-1_c0 <=> Phosphate_c0 + GDP_c0 + H_c0 + Coenzyme_F420_c0');
model = changeGeneAssociation(model,'rxn10525_c0','mmp0937');
%Also add a glutamate to make F420-3, not our mainstream type, but
%something we know that we synthesize (CofF)
model = addReaction(model,{'rxn10527_c0','gamma-F420-2:gamma-L-glutamate ligase'},...
    'L-Glutamate_c0 + GTP_c0 + Coenzyme_F420_c0 <=> Phosphate_c0 + GDP_c0 + H_c0 + Coenzyme_F420-3_c0');
model = changeGeneAssociation(model,'rxn10525_c0','mmp0170');

%Now add the last 2 species to the biomass
[~,f420_idx] = intersect(model.mets,'Coenzyme_F420_c0');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(f420_idx,bio_idx) = -0.0030965;

[~,f420_idx] = intersect(model.mets,'Coenzyme_F420-3_c0');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(f420_idx,bio_idx) = -0.0030965;

%Synthesize oxalate using manual gapfill (see F420 notes)
model = addReaction(model,{'rxn05734_c0','aldehyde dehydrogenase (glyoxylate, NAD)'},...
    'H2O_c0 + NAD_c0 + Glyoxalate_c0 <=> NADH_c0 + 2 H_c0 + Oxalate_c0');

%%%%%%%%%%%%%
%4/6/2015
%%%%%%%%%%%%%
%Remove O2 transporter and exchange; we don't transport O2
model = removeRxns(model,{'rxn05468_c0','EX_cpd00007_e0'});

%%%%%%%%%%%%%
%4/6/2015
%%%%%%%%%%%%%

% Add heme to our list of mets to remove from biomass (above, 2/10/2015)

% Take out the O2 reactions
o2_rxns = findRxnsFromMets(model,'O2_c0');
model = removeRxns(model,o2_rxns);

% Replace the orotate reaction for CTP and UTP synthesis
model = addReaction(model,{'rxn01361_c0','(S)-Dihydroorotate:(acceptor) oxidoreductase'},...
    'NAD_c0 + S-Dihydroorotate_c0 <=> NADH_c0 + H_c0 + Orotate_c0');

%%%%%%%%%%%%%
% 4/14/2015
%%%%%%%%%%%%%
% Model has no way to make fuculose1-phosphate, which is how it gets
% lactate for F420 production. Adding pathway as gapfill:
model = addReaction(model,{'rxn00559_c0','D-Mannose-6-phosphate ketol-isomerase'},...
    'D-mannose-6-phosphate_c0 <=> D-fructose-6-phosphate_c0');
model = addReaction(model,{'rxn00641_c0','GTP:alpha-D-mannose-1-phosphate guanylyltransferase'},...
    'GTP_c0 + D-Mannose1-phosphate_c0 <=> PPi_c0 + GDP-mannose_c0');
model = addReaction(model,{'rxn00642_c0','GDPmannose 4,6-hydro-lyase'},...
    'GDP-mannose_c0 <=> GDP-4-dehydro-D-rhamnose_c0 + H2O_c0');
model = addReaction(model,{'rxn03962_c0','GDP-L-fucose:NADP+ 4-oxidoreductase (3,5-epimerizing)'},...
    'GDP-L-fucose_c0 + NADP_c0 <=> GDP-4-dehydro-D-rhamnose_c0 + NADPH_c0 + H_c0');
model = addReaction(model,{'rxn01431_c0','GTP:L-fucose-1-phosphate guanylyltransferase'},...
    'GTP_c0 + L-Fucose1-phosphate_c0 <=> PPi_c0 + GDP-L-fucose_c0');
model = addReaction(model,{'rxn02262_c0','ATP:6-deoxy-L-galactose 1-phosphotransferase'},...
    'ATP_c0 + L-Fucose_c0 <=> ADP_c0 + L-Fucose1-phosphate_c0');
model = addReaction(model,{'rxn02263_c0','L-Fucose ketol-isomerase'},...
    'L-Fucose_c0 <=> L-Fuculose_c0');
model = addReaction(model,{'rxn02319_c0','ATP:L-fuculose 1-phosphotransferase'},...
    'ATP_c0 + L-Fuculose_c0 <=> ADP_c0 + L-Fuculose1-phosphate_c0');


%%%%%%%%%%%%%
% 4/24/2015
%%%%%%%%%%%%%
% Add EhbN (mmp1153) to the Eha/Ehb reaction:
model = changeGeneAssociation(model,'Eha/Ehb',...
    '(mmp1448 and mmp1449 and mmp1450 and mmp1451 and mmp1452 and mmp1453 and mmp1454 and mmp1455 and mmp1456 and mmp1457 and mmp1458 and mmp1459 and mmp1460 and mmp1461 and mmp1462 and mmp1463 and mmp1464 and mmp1465 and mmp1466 and mmp1467) or (mmp1621 and mmp1622 and mmp1623 and mmp1624 and mmp1625 and mmp1626 and mmp1627 and mmp1628 and mmp1629 and mmp1073 and mmp1074 and mmp1469 and mmp0400 and mmp1153)');

%%%%%%%%%%%%%
% 4/27/2015
%%%%%%%%%%%%%
% Remove the formate-hydrogen lyase, which is misannotated
model = removeRxns(model,'rxn08518_c0');

%%%%%%%%%%%%%
% 4/28/2015
%%%%%%%%%%%%%
% Let Hdr-formate bounds be unrestricted; 3rd H2 -> F420 pathway
model = changeRxnBounds(model,'Hdr_formate',-1000,'l');
model = changeRxnBounds(model,'Hdr_formate',1000,'u');
% Same with HdrABC
model = changeRxnBounds(model,'HdrABC',-1000,'l');
model = changeRxnBounds(model,'HdrABC',1000,'u');
% Same with Eha/Ehb
model = changeRxnBounds(model,'Eha/Ehb',-1000,'l');

% Separate out the CODH and ACS reactions, but keep same genes
% First remove the combined reaction
model = removeRxns(model,'CODH_ACS');
%Add CO-dehydrogenase
model = addReaction(model,{'CODH','Carbon monoxide dehydrogenase'},...
    'CO2_c0 + Reducedferredoxin_c0 + 3 H_c0 <=> CO_c0 + Oxidizedferredoxin_c0 + 0.5 H2_c0 + H2O_c0');
%Associate it with mmpmmp0980,0981,0983,0984,0985
model = changeGeneAssociation(model,'CODH',...
    'mmp0979 and mmp0980 and mmp0981 and mmp0982 and mmp0983 and mmp0984 and mmp0985');

model = addReaction(model,{'ACS','Acetyl-CoA synthase'},...
    '5-Methyl-H4MPT_c0 + CoA_c0 + CO_c0 + 0.5 H2_c0 <=> Acetyl-CoA_c0 + H4MPT_c0 + H_c0');
%Associate it with mmpmmp0980,0981,0983,0984,0985
model = changeGeneAssociation(model,'ACS',...
    'mmp0979 and mmp0980 and mmp0981 and mmp0982 and mmp0983 and mmp0984 and mmp0985');

% Add exchange and diffusion for carbon monoxide
model = addReaction(model,{'rxn10480_c0','CO transporter via diffusion'},...
    'CO_c0 <=> CO_e0');
model = addReaction(model,{'EX_cpd00204_e0','EX_CO_e0'},...
    'CO_e0 <=> ');
% Turn off the supply of CO_e0
model = changeRxnBounds(model,'EX_cpd00204_e0',0,'l');

% Add formulas and charges for F420 things and for CO_c0 and CO_e0
[~,idx] = intersect(model.mets,'CO_c0');
model.metCharge(idx)=1;
model.metFormulas{idx}='CO';
[~,idx] = intersect(model.mets,'CO_e0');
model.metCharge(idx)=1;
model.metFormulas{idx}='CO';
[~,idx] = intersect(model.mets,'Oxalate_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C2O4';
[~,idx] = intersect(model.mets,'7,8-didemethyl-8-hydroxy-5-deazariboflavin_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C16H17N3O7';

[~,idx] = intersect(model.mets,'L-Lactate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C3H5O3';
[~,idx] = intersect(model.mets,'2-phospho-L-lactate_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C3H4O6P';
[~,idx] = intersect(model.mets,'lactyl-(2)-diphospho-(5)-guanosine_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C13H16N5O13P2';
[~,idx] = intersect(model.mets,'Coenzyme_F420-0_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C19H20N3O12P';
[~,idx] = intersect(model.mets,'Coenzyme_F420-1_c0');
model.metCharge(idx)=-3;
model.metFormulas{idx}='C24H26N4O15P';
[~,idx] = intersect(model.mets,'Coenzyme_F420-3_c0');
model.metCharge(idx)=-5;
model.metFormulas{idx}='C34H38N6O21P';
[~,idx] = intersect(model.mets,'GDP-mannose_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C16H23N5O16P2';
[~,idx] = intersect(model.mets,'GDP-4-dehydro-D-rhamnose_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C16H21N5O15P2';
[~,idx] = intersect(model.mets,'GDP-L-fucose_c0');
model.metCharge(idx)=-2;
model.metFormulas{idx}='C16H23N5O15P2';
[~,idx] = intersect(model.mets,'L-Fucose1-phosphate_c0');
model.metCharge(idx)=-1;
model.metFormulas{idx}='C6H12O8P';
[~,idx] = intersect(model.mets,'L-Fucose_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C6H12O5';
[~,idx] = intersect(model.mets,'L-Fuculose_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='C6H12O5';

% Add in F420:NADP oxidoreductase (From genome paper, gene from NCBI,BioCyc/MetaCyc)
model = addReaction(model,{'FNO','F420:NADP oxidoreductase'},...
    'NADP_c0 + Reduced_coenzyme_F420_c0 <=> NADPH_c0 + H_c0 + Coenzyme_F420_c0');
% Associate it with mmpmmp0980,0981,0983,0984,0985
model = changeGeneAssociation(model,'FNO','mmp1550');

%%%%%%%%%%%%%
% 4/30/2015
%%%%%%%%%%%%%
% Restrict reversibility of GAPOR reaction to only be forward (Park 2007)
model = changeRxnBounds(model,'rxn07191_c0',0,'l');

%%%%%%%%%%%%%
% 5/06/2015
%%%%%%%%%%%%%
% Add tranporters

% Add mmp1099 to phosphate transport
model = changeGeneAssociation(model,'rxn05145_c0','(mmp1098 and (mmp0345 or mmp1095) and mmp1097 and mmp1096 and mmp1099) ');

% Add mmp0027 and mmp0484 to potassium transport (TransportDB)
model = changeGeneAssociation(model,'rxn05596_c0','mmp0027 or mmp0484 or mmp1430 or mmp1599');


% Give molybdenum its transporters, then add it to the biomass
model = addReaction(model,{'EX_cpd00131_e0','EX_Molybdenum_e0'},...
    'Molybdenum_e0 <=> ');
model = addReaction(model,{'Mot','ATP-Dependent Molybdenum transport'},...
    'Molybdenum_e0 + H2O_c0 + ATP_c0 <=> Molybdenum_c0 + H_c0 + ADP_c0 + Phosphate_c0');
% Add genes for Mo transport
model = changeGeneAssociation(model,'Mot',...
    '(mmp0205 and mmp0206 and mmp0207) or (mmp0504 and mmp0505 and mmp0506) or (mmp0514 and mmp0515 and mmp0516) or (mmp1650 and mmp1651 and mmp1652)');

% Add Molybdenum to biomass and add its charge/formula
[~,idx] = intersect(model.mets,'Molybdenum_c0');
model.metCharge(idx)=0;
model.metFormulas{idx}='Mo';

[~,idx] = intersect(model.mets,'Molybdenum_e0');
model.metCharge(idx)=0;
model.metFormulas{idx}='Mo';

[~,mo_idx] = intersect(model.mets,'Molybdenum_c0');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(mo_idx,bio_idx) = -0.0030965;

% Add lots of genes to Fe3+ transport (3 gene clusters)
model = changeGeneAssociation(model,'rxn05195_c0',...
    '(mmp0108 and mmp0109 and mmp0110) or (mmp0196 and mmp0197 and mmp0198) or (mmp1181 and mmp1182 and mmp1183)');

% Add mmp1641 for copper transport (TransportDB)
model = changeGeneAssociation(model,'rxn10481_c0',...
    'mmp1165 or mmp1641');

% Add mmp0510 and mmp0720 for calcium transport (TransportDB)
model = changeGeneAssociation(model,'rxn05513_c0',...
    'mmp0510 or mmp0720');

% Add mmp0867 for betaine transport (TransportDB)
model = changeGeneAssociation(model,'rxn05181_c0',...
    'mmp0866 or mmp0867 or mmp0868');

% Add mmp0681 for hypoxanthine and mmp0689 for uracil
model = changeGeneAssociation(model,'rxn05201_c0','mmp0681 or mmp0689');
model = changeGeneAssociation(model,'rxn05197_c0','mmp0681 or mmp0689');

% Add mmp0867 for proline transport (TransportDB)
model = changeGeneAssociation(model,'rxn05165_c0',...
    'mmp0866 or mmp0867 or mmp0868');

% Add mmp1197 and mmp1198 for nitrate transport and taurine transport (TransportDB)
model = changeGeneAssociation(model,'rxn05171_c0',...
    'mmp1260 or mmp1197 or mmp1198');
model = changeGeneAssociation(model,'rxn05172_c0',...
    'mmp1260 or mmp1197 or mmp1198');

% Alter genes for Na+/H+ antiporter
model = changeGeneAssociation(model,'rxn05209_c0',...
    'mmp0100 or mmp0587 or mmp0679 or mmp0707 or mmp0864');

%%%%%%%%%%%%%
% 6/04/2015
%%%%%%%%%%%%%
% Add an ID to the model
model.id = 'M_maripaludis_S2_genome_scale';

%%%%%%%%%%%%%
% 7/08/2015
%%%%%%%%%%%%%
% Add homocysteine synthesis reaction from literature
model = addReaction(model,{'HcyS','Homocysteine Synthesis'},...
    'L-Aspartate4-semialdehyde_c0 + H2S_c0 + H2_c0 <=> Homocysteine_c0 + H2O_c0');
model = changeGeneAssociation(model,'HcyS','mmp1358 and mmp1359');

% Remove reactions based on the 2015 Allen paper on homocysteine
%model = removeRxns(model,{'rxn01304_c0','rxn05957_c0'});

% Revise last step of methionine synthesis (Remove folate steps)
model = removeRxns(model,...
    {'rxn01304_c0','rxn03004_c0','rxn03137_c0','rxn00692_c0','rxn04954_c0','rxn11944_c0'});

% Add the step for methionine synthesis
model = addReaction(model,{'MetS','Methionine Synthase'},...
    'Homocysteine_c0 + 5-Methyl-H4MPT_c0 -> L-Methionine_c0 + H4MPT_c0');

% Add subsystems 
[~, idx] = intersect(model.rxns,'HcyS');
model.subSystems{idx} = 'Methionine Biosynthesis';
[~, idx] = intersect(model.rxns,'MetS');
model.subSystems{idx} = 'Methionine Biosynthesis';

% Make SAICAR reaction reversible to fix Histidine production
model = changeRxnBounds(model,'rxn03147_c0',-1000,'l');

% Try to fix the folates missing by replacing their reactions with
% corresponding methanopterins
% We need this purely to consume the Dihydrofolate
%rxn00686_c0	H_c0 + NADPH_c0 + Dihydrofolate_c0 	->	NADP_c0 + Tetrahydrofolate_c0 	
model = addReaction(model,{'rxn02430_c0','Dihydromethanopterin reductase'},...
    'H_c0 + NADPH_c0 + 7,8-Dihydromethanopterin_c0 <=> NADP_c0 + H4MPT_c0');

% Add this too for the hydropantoate, which is needed for CoA synthesis
%rxn00912_c0	H2O_c0 + 3-Methyl-2-oxobutanoate_c0 + 5-10-Methylenetetrahydrofolate_c0 	->	2-Dehydropantoate_c0 + Tetrahydrofolate_c0
model = addReaction(model,{'H4MPT3M2Om','5,10-Methylenetetrahydromethanopterin 3-Methyl-2-oxobutanoate methyltransferase'},...
    'H2O_c0 + 3-Methyl-2-oxobutanoate_c0 + 5_10-Methylenetetrahydromethanopterin_c0 <=>	2-Dehydropantoate_c0 + H4MPT_c0');

% Might not need this either, but it's probably important for dUMP and
% dTMP, so add it and leave the Dihydrofolate for now
%rxn01520_c0	dUMP_c0 + 5-10-Methylenetetrahydrofolate_c0 	->	dTMP_c0 + Dihydrofolate_c0 	(mmp1379 or mmp0986)
model = addReaction(model,{'H4MPTdUMPm','5,10-Methylenetetrahydromethanopterin dUMP C-methyltransferase'},...
    'dUMP_c0 + 5_10-Methylenetetrahydromethanopterin_c0	<=> dTMP_c0 + 7,8-Dihydromethanopterin_c0');

%rxn03137_c0	AICAR_c0 + 10-Formyltetrahydrofolate_c0 	->	FAICAR_c0 + Tetrahydrofolate_c0 	
model = addReaction(model,{'FH4MPTAf','Formyl-H4MPT AICAR Formyltransferase'},...
    'AICAR_c0 + 5-Formyl-H4MPT_c0 <=> FAICAR_c0 + H4MPT_c0');

% Remove reaction with Methylene-tetrahydrofolate or Dihydrofolate
rxns = findRxnsFromMets(model,{'5-10-Methylenetetrahydrofolate_c0','Dihydrofolate_c0'});
model = removeRxns(model,rxns);

% Add info for Dihydromethanopterin
[~,idx] = intersect(model.mets,'7,8-Dihydromethanopterin_c0');
model.metCharge(idx) = -3;
model.metFormulas{idx} = 'C30H40N6O16P';
model.metSEEDID{idx} = 'cpd03523';

%%%%%%%%%%%%%
% 7/14/2015
%%%%%%%%%%%%%
% Remove "SINK" and replace it with a hypothetical hydroxypyruvaldehyde
% reduction to G3P
model = removeRxns(model,'SINK');
model = addReaction(model,{'HPAr','Hydroxypyruvaldehyde phosphate reductase'},...
    'Hydroxypyruvaldehyde_phosphate_c0 + NADH_c0 <=> Glyceraldehyde3-phosphate_c0 + NAD_c0');

%Also, remove reactions that have folate
rxns = findRxnsFromMets(model,{'Folate_c0','Folate_e0'});
model = removeRxns(model,rxns);


%%%%%%%%%%%%%
%7/15/2015
%%%%%%%%%%%%%

% Add final step of acetamido sugar things (*Still need a consumption step)
model = addReaction(model,{'rxn02377_c0','UDP-N-acetyl-D-mannosamine:NAD+ 1-oxidoreductase'},...
    'H2O_c0 + 2 NAD_c0 + UDP-N-acetyl-D-mannosamine_c0 <=> 2 NADH_c0 + 3 H_c0 + UDP-N-acetyl-D-mannosaminouronate_c0');
model = changeGeneAssociation(model,'rxn02377_c0','mmp0706');

% Add info for the metabolite we just added
[~,idx] = intersect(model.mets,'UDP-N-acetyl-D-mannosaminouronate_c0');
model.metCharge(idx) = -3;
model.metFormulas{idx} = 'C17H22N3O18P2';
model.metSEEDID{idx} = 'cpd03732';

% Add subsystem info for reaction
[~, idx] = intersect(model.rxns,'rxn02377_c0');
model.subSystems{idx} = 'None';

%%%%%%%%%%%%%
% 4/16/2015
%%%%%%%%%%%%%
% Remove dead ends that have no genes
model = removeDeadGapFills(model);

%%%%%%%%%%%%%
%9/19/2014
%%%%%%%%%%%%%
%One of final steps should always be to add the kbase aliases:
model = addKbaseAliases(model);

%%%%%%%%%%%%%
%7/14/2015
%%%%%%%%%%%%%

% Last step: fix some reaction names that aren't right and add transport
% for things that need it

% Change name of CH4 from 'EX_Methane_c0' to 'Methane_c0' and 'Methane_e0'
[~,idx] = intersect(model.mets,'EX_Methane_c0');
model.mets{idx} = 'Methane_c0';
model.metNames{idx} = 'Methane_c0';

% Add a transport (diffusion) for both hydrogen and methane
model = addReaction(model,{'rxn10542_c0','Hydrogen transport via diffusion'},...
    'H2_e0 <=> H2_c0');
[~,idx] = intersect(model.rxns,'rxn10542_c0');
model.subSystems{idx} = 'Transport';

model = addReaction(model,{'rxn10471_c0','Methane transport via diffusion'},...
    'Methane_e0 <=> Methane_c0');
[~,idx] = intersect(model.rxns,'rxn10471_c0');
model.subSystems{idx} = 'Transport';

% Fix the names and formulas of their exchanges
% Methane
[~,idx] = intersect(model.rxns,'Ex_cpd01024_c0');
model.rxns{idx} = 'EX_cpd01024_e0';
model = addReaction(model,{'EX_cpd01024_e0','EX_Methane_e0'},...
    'Methane_e0 <=> ');

[~,idx] = intersect(model.rxns,'Ex_cpd11640_c0');
model.rxns{idx} = 'EX_cpd11640_e0';
model = addReaction(model,{'EX_cpd11640_e0','EX_Hydrogen_e0'},...
    'H2_e0 <=> ');

% Fix names of Cob(I)yrinate_diamide and Cob(II)yrinate_diamide
[~,idx ] = intersect(model.mets,'CobIyrinate_diamide_c0');
model.mets{idx} = 'Cob(I)yrinate_diamide_c0';
model.metNames{idx} = 'Cob(I)yrinate_diamide_c0';
[~,idx ] = intersect(model.mets,'CobIIyrinate_diamide_c0');
model.mets{idx} = 'Cob(II)yrinate_diamide_c0';
model.metNames{idx} = 'Cob(II)yrinate_diamide_c0';

%%%%%%%%%%%%%
%7/16/2015
%%%%%%%%%%%%%
% Very last step: add free energy values for 1 mM from Equilibrator site
model = addDG2MM(model);
