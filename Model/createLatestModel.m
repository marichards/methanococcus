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
%5/02/2014 (All changes but changing rxnNames)
%5/22/2014 (Change IDs to abbreviations, add rxnNames)
%%%%%%%%%%%%%%%%%%
%Add the 25 reactions additionally marked for addition
%2 Reactions from "other modified reactions"
%Modify the charge balance for the first one
model = addReaction(model,'rxn06874_c0',...
    '16 H2O_c0 + 16 ATP_c0 + 2 H_c0 + 8 Reducedferredoxin_c0 + N2_c0 <=> 16 Phosphate_c0 + 16 ADP_c0 + 2 NH3_c0 + 8 Oxidizedferredoxin_c0 + H2_c0');
model = addReaction(model,'rxn04042_c0',...
    'ADP_c0 + D-Glucose_c0 <=> D-glucose-6-phosphate_c0 + AMP_c0');
model.rxnNames{end}='ADP:D-glucose 6-phosphotransferase';
model = addReaction(model,'rxn04043_c0',...
    'ADP_c0 + D-fructose-6-phosphate_c0 <=> D-fructose-1_6-bisphosphate_c0 + AMP_c0');
model.rxnNames{end}='ADP:D-fructose-6-phosphate 1-phosphotransferase';


%21 Reactions from "added reactions"
model = addReaction(model,'rxn00340_c0',...
    'ATP_c0 + L-Aspartate_c0 + NH3_c0 -> PPi_c0 + AMP_c0 + L-Asparagine_c0');
model.rxnNames{end}='L-Aspartate:ammonia ligase (AMP-forming)';
model = addReaction(model,'rxn00184_c0',...
    '2-Oxoglutarate_c0 + H_c0 + NADPH_c0 + NH3_c0 <=> L-Glutamate_c0 + H2O_c0 + NADP_c0');
model.rxnNames{end}='L-Glutamate:NADP+ oxidoreductase (deaminating)';
model = addReaction(model,'rxn10603_c0',...
    'L-Cysteine_c0 + sulfoacetaldehyde_c0 <=> H2O_c0 + 2-(sulfomethyl)thiazolidine-4-carboxylate_c0');
model.rxnNames{end}='Sulfomethyl thiazolidine synthase';
model = addReaction(model,'rxn10598_c0',...
    '2-(sulfomethyl)thiazolidine-4-carboxylate_c0 + NADH_c0 + H_c0 <=> sulfoethylcysteine_c0 + NAD_c0');
model.rxnNames{end}='Sulfoethylcysteine synthase';
model = addReaction(model,'rxn10479_c0',...
    'H2O_c0 + sulfoethylcysteine_c0 <=> NH3_c0 + CoM_c0 + Pyruvate_c0');
model.rxnNames{end}='Coenzyme M synthase';
model = addReaction(model,'rxn00974_c0',...
    'Citrate_c0  <=> H2O_c0 + cis-Aconitate_c0');
model.rxnNames{end} ='Citrate hydrolase';
model = addReaction(model,'rxn01388_c0',...
    'H2O_c0 + cis-Aconitate_c0  <=> D-threo-Isocitrate_c0');
model.rxnNames{end} ='Isocitrate hydrolase';
model = addReaction(model,'rxn00198_c0',...
    'NADP_c0 + D-threo-Isocitrate_c0  <=> CO2_c0 + NADPH_c0 + 2-oxoglutarate_c0');
model.rxnNames{end} ='Isocitrate:NADP+ oxidoreductase (decarboxylating)';
model = addReaction(model,'rxn08518_c0',...
    'Formate_c0 + H_c0  -> CO2_c0 + H2_c0');
model.rxnNames{end} ='Formate-hydrogen lyase';
model = addReaction(model,'rxn03079_c0',...
    'H_c0 + Coenzyme_F420_c0 + 5_10-Methylenetetrahydromethanopterin_c0  <=> H2_c0 + Reduced_coenzyme_F420_c0 + 5_10-Methenyltetrahydromethanopterin_c0');
model = addReaction(model,'rxn00256_c0',...
    'Citrate_c0 + CoA_c0 + H_c0 <=> Acetyl-CoA_c0 + H2O_c0 + Oxaloacetate_c0');
model.rxnNames{end} ='Citrate oxaloacetate-lyase';
model = addReaction(model,'rxn04026_c0',...
    '3-sulfopyruvate_c0 + H_c0 -> CO2_c0 + sulfoacetaldehyde_c0');
model.rxnNames{end} ='Sulfopyruvate carboxy-lyase';
%Modified (R)-sulfolactate to (2R)-3-sulfolactate to match the seed on
%5/20/2014
model = addReaction(model,'rxn04934_c0',...
    '(2R)-3-sulfolactate_c0 + NAD_c0 <=> NADH_c0 + H_c0 + 3-sulfopyruvate_c0');
model.rxnNames{end} ='(R)-2-hydroxyacid:NAD+ oxidoreductase';
%Modified (R)-sulfolactate to (2R)-3-sulfolactate to match the seed on
%5/20/2014
model = addReaction(model,'rxn04036_c0',...
    '2R-Phosphosulfolactate_c0 + H2O_c0 -> Phosphate_c0 + H_c0 + (2R)-3-sulfolactate_c0');
model.rxnNames{end} ='(R)-2-phospho-3-sulfolactate phosphohydrolase';
model = addReaction(model,'rxn07741_c0',...
    'tetrahydrodipicolinate_c0 + H_c0 + H2O_c0 + L-Glutamate_c0 <=> 2-oxoglutarate_c0 + LL-2_6-Diaminopimelate_c0');
model.rxnNames{end} ='LL-2,6-diaminoheptanedioate:2-oxoglutarate aminotransferase';
model = addReaction(model,'rxn05109_c0',...
    'Acetyl-CoA_c0 + H2O_c0 + Pyruvate_c0 -> H_c0 + CoA_c0 + Citramalate_c0');
model.rxnNames{end} ='Citramalate_synthase';
model = addReaction(model,'rxn02749_c0',...
    'Citramalate_c0 <=> H2O_c0 + Citraconate_c0');
model.rxnNames{end} ='(R)-2-Methylmalate hydro-lyase';
%In following 2 reactions, replaced beta-methyl-d-malate with D-erythro-3-methylmalate
model = addReaction(model,'rxn02751_c0',...
    'H2O_c0 + Citraconate_c0 <=> D-erythro-3-methylmalate_c0');
model.rxnNames{end} ='Isopropylmalate_isomerase';
model = addReaction(model,'rxn00735_c0',...
    'D-erythro-3-methylmalate_c0 + NAD_c0 <=> NADH_c0 + CO2_c0 + 2-Oxobutyrate_c0');
model.rxnNames{end} ='Methylmalate_dehydrogenase';
model = addReaction(model,'rxn08043_c0',...
    '2-Oxobutyrate_c0 + H_c0 + Pyruvate_c0 -> 2-Aceto-2-hydroxybutanoate_c0 + CO2_c0');
model.rxnNames{end} ='2-aceto-2-hydroxybutanoate synthase';
model = addReaction(model,'rxn08764_c0',...
    'NADPH_c0 + H_c0 + 2-Aceto-2-hydroxybutanoate_c0 <=> 2_3-Dihydroxy-3-methylvalerate_c0 + NADP_c0');
model.rxnNames{end} ='ketol-acid reductoisomerase (2-Acetolactate)';

%Associate genes with added reactions
model = changeGeneAssociation(model,'rxn06874_c0','mmp0853 or mmp0856 or mmp0857');
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
model = addReaction(model,'rxn00405_c0',...
    'L-Arginine_c0 + H_c0 -> CO2_c0 + Agmatine_c0 ');
model.rxnNames{end} = 'Arginine carboxy-lyase';
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
model = addReaction(model,'rxn00171_c0',...
    'Acetaldehyde_c0 + NAD_c0 + CoA_c0 <=> H_c0 + NADH_c0 + Acetyl-CoA_c0 ');
model.rxnNames{end} = 'Acetaldehyde:NAD+ oxidoreductase (CoA-acetylating';
%Associate the gene 'mmp1423' with reaction 'aldehyde dehydrogenase'
model = changeGeneAssociation(model,'rxn00171_c0','mmp1423');

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
model = addReaction(model,'CODH',...
    'CO2_c0 + CoA_c0 + 2 H_c0 + Reducedferredoxin_c0 + 5-Methyl-H4MPT_c0 <=> Acetyl-CoA_c0 + Oxidizedferredoxin_c0 + H2O_c0 + H4MPT_c0');
model.rxnNames{end} ='CO dehydrogenase';
%Associate it with mmpmmp0980,0981,0983,0984,0985
model = changeGeneAssociation(model,'CODH',...
    'mmp0979 and mmp0980 and mmp0981 and mmp0982 and mmp0983 and mmp0984 and mmp0985');
%Add acetate exchange and transport, both reversible
model = addReaction(model,'EX_cpd00029_e0',...
    'Acetate_e0 <=> ');
model.rxnNames{end} = 'Acetate exchange';
model = addReaction(model,'ACT',...
    'Acetate_e0 <=> Acetate_c0');
model.rxnNames{end} = 'Acetate transport';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of 5/15 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5/20 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add charges for added metabolites
%There are 13 added metabolites...add charges 1 by 1 (from Kbase)
%    'sulfoacetaldehyde_c0'
model.metCharge(end-12)=-1;
%    '2-(sulfomethyl)thiazolidine-4-carboxylate_c0'
model.metCharge(end-11)=-1;
%    'sulfoethylcysteine_c0'
model.metCharge(end-10)=-1;
%    'Citrate_c0'
model.metCharge(end-9)=-3;
%    'cis-Aconitate_c0'
model.metCharge(end-8)=-3;
%    'D-threo-Isocitrate_c0'
model.metCharge(end-7)=-3;
%    '2-oxoglutarate_c0'
model.metCharge(end-6)=-2;
%    '3-sulfopyruvate_c0'
model.metCharge(end-5)=-2;
%    '(2R)-3-sulfolactate_c0'
model.metCharge(end-4)=-2;
%    'Citramalate_c0'
model.metCharge(end-3)=-2;
%    'Citraconate_c0'
model.metCharge(end-2)=-2;
%    'D-erythro-3-Methylmalate_c0'
model.metCharge(end-1)=-2;
%    'Acetate_e0'
model.metCharge(end)=-1;

%Fix charges for 3 other reactions
%Reaction 07191_c0; change 2 Fd to 1
model = addReaction(model,'rxn07191_c0',...
    'H2O_c0 + Glyceraldehyde3-phosphate_c0 + Oxidizedferredoxin_c0 <=> 3.000000 H_c0 + 3-Phosphoglycerate_c0 + Reducedferredoxin_c0');
%Reaction 04045_c0; Remove 2 protons from the left
model = addReaction(model,'rxn04045_c0',...
    'Sirohydrochlorin_c0 + Co2_c0 	<=>	Cobalt-precorrin_2_c0');
%Reaction 05029_c0; add 2 protons to the right
model = addReaction(model,'rxn05029_c0',...
    'ATP_c0 + Cobinamide_c0 ->	Triphosphate_c0 + Adenosyl_cobinamide_c0 + 2 H_c0');

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

%change the gene association of rxn00371_c0 from mmp1233 to (mmp0138 or mmp0139 or mmp1297 or mmp1298)
model = changeGeneAssociation(model,'rxn00371_c0','mmp0138 or mmp0139 or mmp1297 or mmp1298');

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
%EX Urea e0	-1000
model = removeRxns(model,'EX_cpd00073_e0');
%EX fe3 e0	-1000
model = removeRxns(model,'EX_cpd10516_e0');
%EX Spermine e0	-1000
model = removeRxns(model,'EX_cpd00558_e0');
%EX Nitrate e0	-1000
model = removeRxns(model,'EX_cpd00209_e0');
%EX BET e0	-1000
model = removeRxns(model,'EX_cpd00540_e0');
%EX Cytosine e0	-1000
model = removeRxns(model,'EX_cpd00307_e0');
%EX glycogenn-1 c0	-1000
model = removeRxns(model,'EX_cpd15302_c0');
%EX Uracil e0	-1000
model = removeRxns(model,'EX_cpd00092_e0');
%EX ddca e0	-1000
model = removeRxns(model,'EX_cpd01741_e0');
%EX Dephospho-CoA e0	-1000
model = removeRxns(model,'EX_cpd00655_e0');
%EX Cobinamide e0	-1000 %Necessary for Colamide?
model = removeRxns(model,'EX_cpd03422_e0');
%EX Oxidized glutathione e0	-1000
model = removeRxns(model,'EX_cpd00111_e0');

%Metals in the biomass that must be removed
%For now, don't remove their transporters or exchanges
%EX Zn2 e0	-1000
model = changeRxnBounds(model,'EX_cpd00034_e0',0,'l');
%EX Cu2 e0	-1000
model = changeRxnBounds(model,'EX_cpd00058_e0',0,'l');
%EX Mn2 e0	-1000
model = changeRxnBounds(model,'EX_cpd00030_e0',0,'l');

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
model = addReaction(model,'rxn00781_c0',...
    'Phosphate_c0 + NADP_c0 + Glyceraldehyde3-phosphate_c0 <=> NADPH_c0 + 1_3-Bisphospho-D-glycerate_c0');

%Change Gene Associations to AND for ATP Synthase and rxn00250_c0, plus add
%some genes to ATP Synthase
model = changeGeneAssociation(model,'rxn08173_c0',...
    'mmp1038 and mmp1039 and mmp1040 and mmp1041 and mmp1042 and mmp1043 and mmp1044 and mmp1045 and mmp1046');
model = changeGeneAssociation(model,'rxn00250_c0','mmp0340 and mmp0341');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6/23 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Turn off  rxn07191_c0, make rxn05398_c0 and CODH irrev (John spreadsheet 6/19)
%model = changeRxnBounds(model,'rxn07191_c0',0,'b');
%model = changeRxnBounds(model,'rxn05938_c0',0,'l');
%model = changeRxnBounds(model,'CODH',0,'l');

%Add genes for Eha: mmp1448-1467
model = changeGeneAssociation(model,'Eha',...
    'mmp1448 and mmp1449 and mmp1450 and mmp1451 and mmp1452 and mmp1453 and mmp1454 and mmp1455 and mmp1456 and mmp1457 and mmp1458 and mmp1459 and mmp1460 and mmp1461 and mmp1462 and mmp1463 and mmp1464 and mmp1465 and mmp1466 and mmp1467');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6/26 model changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Group EhA, HdrABC, and rxn11938_c0 together by giving them a specific
%ferredoxin
model = addReaction(model,'Eha',...
    'Fdox*1_c0 + 2.000000 Na_e0 + H2_c0 ->	2.000000 H_c0 + Fdred*1_c0 + 2.000000 Na_c0');
model = addReaction(model,'HdrABC',...
    'Fdox*1_c0 + CoM-S-S-CoB_c0 + 2.000000 H2_c0 ->	2.000000 H_c0 + Fdred*1_c0 + CoM_c0 + HTP_c0'); 	
model = addReaction(model,'rxn11938_c0',...
    'H2O_c0 + Fdox*1_c0 + Formylmethanofuran_c0	<=>	H_c0 + CO2_c0 + Fdred*1_c0 + Methanofuran_c0');

%Add EhB, indolepyruvate oxidoreductase with new specific ferredoxin; both are dead for now
model = addReaction(model,'Ehb',...
    'Fdox*2_c0 + 2.000000 Na_e0 + H2_c0 ->	2.000000 H_c0 + Fdred*2_c0 + 2.000000 Na_c0');
model = addReaction(model,'IPOR',...
    'Indole-3-pyruvate_c0 + Fdox*2_c0 + CoA_c0 <=> S-2-(indol-3 yl)acetyl-CoA_c0 + CO2_c0 + Fdred*2_c0');
%Give the CODH, rxn05938, and rxn05939 the same ferredoxin
model = addReaction(model,'CODH',...
    'CO2_c0 + CoA_c0 + 2 H_c0 + Fdred*2_c0 + 5-Methyl-H4MPT_c0 <=> Acetyl-CoA_c0 + Fdox*2_c0 + H2O_c0 + H4MPT_c0');
model = addReaction(model,'rxn05938_c0',...
    'H_c0 + CO2_c0 + Acetyl-CoA_c0 + Fdred*2_c0 <=>	CoA_c0 + Fdox*2_c0 + Pyruvate_c0');
model = addReaction(model,'rxn05939_c0',...
    'H_c0 + CO2_c0 + Fdred*2_c0 + Succinyl-CoA_c0 	<=>	CoA_c0 + 2-Oxoglutarate_c0 + Fdox*2_c0');

%Add Genes for HdrABC and take them from rxn03126_c0
model = changeGeneAssociation(model,'HdrABC',...
    '(mmp0825 or mmp1607) and ((mmp1053 and mmp1054) or (mmp1155 and mmp1154))');
model = changeGeneAssociation(model,'rxn03126_c0','');

%Add Genes for EhB and indolepyruvate oxidoreductase
model = changeGeneAssociation(model,'Ehb',...
    'mmp1621 and mmp1622 and mmp1623 and mmp1624 and mmp1625 and mmp1626 and mmp1627 and mmp1628 and mmp1629 and mmp1073 and mmp1074 and mmp1469 and mmp400');
model = changeGeneAssociation(model,'IPOR',...
    '(mmp315 and mmp316) or (mmp713 and mmp714)');