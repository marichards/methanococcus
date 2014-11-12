function [recon_rxns,gapfill_rxns]=createRxnConfidenceSheet(model,workbook,sheet)
%
%Create an Excel sheet that includes whether each reaction was in the
%original reconstruction, gap-filled to make the original model, or
%manually added later on 


%Load the original model (M_mar)
load('original_model.mat')

%Remove the non-genes from the original model
M_mar=removeGene(M_mar,'Unknown');
M_mar=removeGene(M_mar,'fig');

%Find the reactions that had genes in the orignal model
recon_rxns = setdiff(M_mar.rxns,findRxnsWOGenes(M_mar));

%Find transport and exchange reactions, give them their own tags
exc_rxns = M_mar.rxns(findExcRxns(M_mar));

%Assign the handfull of reactions we added so they're pulled out
manual_rxns = {'ATP_synthase';'HdrABC';'Eha'};

%Find the rest, assign them to the "GapFilled" thing
gapfill_rxns = setdiff(setdiff(setdiff(M_mar.rxns,recon_rxns),exc_rxns),manual_rxns);

%We now know the origin of all things in the original. Now check which of
%those are actually in the new reaction and add that tag accordingly

%Step 1: Create a set of tags call "tags"
tags = cell(length(model.rxns),1);

%Step 2: Iterate through each of the first 3 lists, find the reactions in the model, and
%then add the appropriate tag for each

%Use 'OR' for 'Original Reconstruction'
[~,idx] = intersect(model.rxns,recon_rxns);
for i=1:length(idx)
    tags{idx(i)}='OR';
end

%Use 'EX' for 'Exchange Reaction'
[~,idx] = intersect(model.rxns,exc_rxns);
for i=1:length(idx)
    tags{idx(i)}='EX';
end

%Use 'GF' for 'Gapfilled Reaction'
[~,idx] = intersect(model.rxns,gapfill_rxns);
for i=1:length(idx)
    tags{idx(i)}='GF';
end

%Step 3: Basically add "AM" for reactions that were "Added Manually"
[~,idx] = setdiff(model.rxns,M_mar.rxns);
for i=1:length(idx)
    tags{idx(i)}='AM';
end

%Also do it for anything that's in the manual set
[~,idx] = intersect(model.rxns,manual_rxns);
for i=1:length(idx)
    tags{idx(i)}='AM';
end

%Print out the tags and other reaction info into the excel sheet designated
formulas = printRxnFormula(model,model.rxns,'False');
A = [{'RxnId','RxnName','Origin Tag','Gene-Rxn Rule','Rxn Formula'};...
    model.rxns,model.rxnNames,tags,model.grRules,formulas];
xlswrite(workbook,A,sheet);
