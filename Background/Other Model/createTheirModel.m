function theirs = createTheirModel()

%The model file we grabbed is faulty, so I'm going to recreate the model
%using the excel file 'supplementary 1.xls'

%Reactions are in sheet "Reactions" with rxns in A2:572, names in B2:572,
%formulas in C2:572, GPRs in D2:564, Proteins in E3:562, Subsystems in
%F2:571, and ECs in H2:570
[~,rxnAbrList] = xlsread('supplementary 1.xls','Reactions','A2:A572');
[~,rxnNameList] = xlsread('supplementary 1.xls','Reactions','B2:B572');

%Exchange Reactions are in the sheet "exchange reactions" with names in
%A2:49 and formulas in B2:49

%Metabolites should be added automatically, but perhaps we can just add
%them instead....? Sheet "Metabolites", all 2:606, A is met, B is name, C
%is formula, E is Kegg ID

%Add minimal media from "Minimal media" sheet

 model = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList,...
    lowerBoundList,upperBoundList,subSystemList,grRuleList,geneNameList,...
    systNameList)

INPUTS
 rxnAbrList            List of names of the new reactions
 rxnNameList           List of names of the new reactions
 rxnList               List of reactions: format: {'A -> B + 2 C'}
                       If the compartment of a metabolite is not
                       specified, it is assumed to be cytoplasmic, i.e. [c]

OPTIONAL INPUTS
 revFlagList           List of reversibility flag (opt, default = 1)
 lowerBoundList        List of lower bound (Default = 0 or -vMax)
 upperBoundList        List of upper bound (Default = vMax)
 subSystemList         List of subsystem (Default = '')
 grRuleList            List of gene-reaction rule in boolean format (and/or allowed)
                       (Default = '');
 geneNameList          List of gene names (used only for translation
                       from common gene names to systematic gene names)
 systNameList          List of systematic names

OUTPUT
 model                 COBRA model structure