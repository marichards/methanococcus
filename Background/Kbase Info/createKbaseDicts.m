function [compounds,reactions] = createKbaseDicts()

%Using the excel sheets, we want to pull out info...

%Pull out the compound IDs and KEGG compound IDs from the compounds along
%with the names
%Goes from row 2-16280
%SEED IDs are in column A, names are column B, KEGG IDs are column E
%First bring in SEED IDs
[~,txt] = xlsread('ModelSEED-compounds-db.xls','Compounds','A2:A16280');
%These are full, so just put them in
compounds.seedIDs = txt;

%Then bring in Names
[~,txt] = xlsread('ModelSEED-compounds-db.xls','Compounds','B2:B16280');
%These are also full, so bring them in
compounds.names = txt;

%Then bring in KEGG IDs
[~,txt] = xlsread('ModelSEED-compounds-db.xls','Compounds','E2:E16280');
%These are missing 5 entries. Add them
txt = [txt;{'';'';'';'';''}];
compounds.keggIDs = txt;

%Then bring in free energy of formation (column H)
%[~,txt] = xlsread('ModelSEED-compounds-db.xls','Compounds','H2:H16280');
%compounds.dGs = txt;

%Pull out the KEGG Rxn IDs, EC numbers...subsystems aren't there
%Goes from row 2-13273
%Rxn IDs are A, ECs are C, Kegg are D, name are B (don't need but take)

%Bring in Seed IDs
[~,txt] = xlsread('ModelSEED-reactions-db.xls','Reactions','A2:A13273');
%These are full, bring them in
reactions.seedIDs = txt;

%Bring in Seed Names
[~,txt] = xlsread('ModelSEED-reactions-db.xls','Reactions','B2:B13273');
%These are full, bring them in
reactions.names = txt;

%Bring in EC numbers
[~,txt] = xlsread('ModelSEED-reactions-db.xls','Reactions','C2:C13273');
%These are missing 6, 3 at the top, 3 at the bottom
txt = [{'';'';''};txt;{'';'';''}];
reactions.ECs = txt;

%Bring in Kegg IDs
[~,txt] = xlsread('ModelSEED-reactions-db.xls','Reactions','D2:D13273');
%These are missing 17, 1 at the top, 16 at the bottom
txt = [{''};txt;{'';'';'';'';'';'';'';'';'';'';'';'';'';'';'';''}];
reactions.keggIDs = txt;
