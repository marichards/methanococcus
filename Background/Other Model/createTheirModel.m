function theirs = createTheirModel()

%The model file we grabbed is faulty, so I'm going to recreate the model
%using the excel file 'supplementary 1.xls'

%Reactions are in sheet "Reactions" with rxns in A2:572, names in B2:572,
%formulas in C2:572, GPRs in D2:564, Proteins in E3:562, Subsystems in
%F2:571, and ECs in H2:570
[~,rxnAbrList] = xlsread('supplementary 1.xls','Reactions','A2:A572');
[~,rxnNameList] = xlsread('supplementary 1.xls','Reactions','B2:B572');
[~,rxnList] = xlsread('supplementary 1.xls','Reactions','C2:C572');

[~,excRxnList] = xlsread('supplementary 1.xls','exchange reactions',...
    'A2:A49');
excNameList = excRxnList;
[~,excList] = xlsread('supplementary 1.xls','exchange reactions','B2:B49');

%Put them together
rxnAbrList = [rxnAbrList;excRxnList];
rxnNameList = [rxnNameList;excNameList];
rxnList = [rxnList;excList];

%Grab the reactions
grRuleList = xlsread('supplementary 1.xls','Reactions','D2:D564');
grRuleList = [grRuleList;{'';'';'';'';'';'';'';''}];

%Exchange Reactions are in the sheet "exchange reactions" with names in
%A2:49 and formulas in B2:49

%Metabolites should be added automatically, but perhaps we can just add
%them instead....? Sheet "Metabolites", all 2:606, A is met, B is name, C
%is formula, E is Kegg ID

%Create the model using what we read in
theirs = createModel(rxnAbrList,rxnNameList,rxnList);%,[],[],[],{},grRuleList);

%Add GPRs
[~,grRuleList] = xlsread('supplementary 1.xls','Reactions','D2:D564');
grRuleList = [grRuleList];
for i = 1:length(grRuleList)
    
    %Add the gene association rule and make it lower case
    theirs = changeGeneAssociation(theirs,theirs.rxns{i},lower(grRuleList{i}));
end

%Add other info not in the model
[~,ECs] = xlsread('supplementary 1.xls','Reactions','H2:H572');
theirs.rxnECNumbers = [ECs;cell(50,1)];

%[~,metKEGGID] = xlsread('supplementary 1.xls','Metabolites','E2:E606');

%Add bounds
[theirs,unchanged] = addBounds(theirs);


%Exchanges should be 0 lb except for media (no hydron (H+), leave out)
media = {'Hydrogen exchange';...
        'Hydrogen sulfide exchange';...
        'ammonium ion exchange';...
        'Orthophosphate exchange';...
        'sulfite exchange';...
        'Carbon dioxide exchange';...
        'proton exchange';...
        'ferrous ion exchange';...
        'water exchange';...
        'Nickel exchange';...
        'cobalt exchange';...
        'Tungstate exchange';...
        'sodium exchange'...
        };
    
 %Except these from the bounds we want to make 0 on lower
 %change_lower = setdiff(unchanged,media);
 
 %Change the lower bounds
 %theirs = changeRxnBounds(theirs,change_lower,0,'l');
 
 %Add their objective (biomass)
 theirs = changeObjective(theirs,'R571',1);
 
 %Remove the other biomass exchange, it's not a thing
 theirs = removeRxns(theirs,'Biomass Formation');
 
        

        