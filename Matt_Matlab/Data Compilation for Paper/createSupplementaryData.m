function createSupplementaryData(model)

% Takes in the M. maripaludis model and creates the automated portion of
% the Supplementary Materials Excel workbook containing: 
% Reactions: info for all reactions, including cross references, formulas,
% GPRs, subsystems, names, bounds,free energy (if applicable)
% Metabolites: info for all metabolites, including cross references,
% formulas, charges
% FVA Results on H2: bounds for reaction fluxes on H2+CO2 media
% FVA Results on Formate: bounds for reaction fluxes on formate media
% 

% Compile the reactions information
formulas = printRxnFormula(model,model.rxns,false,false,true);
% Add the confidence data from other code (below)
tags = createRxnConfidenceSheet(model);
A = [{'ID','Name','KEGG ID','Formula','EC Number(s)','GPR Rules',...
    'Subsystem','Lower Bound','Upper Bound','Free Energy (kJ/mmol)',...
    'Origin Tag'};...
    model.rxns,model.rxnNames,model.rxnKEGGID,formulas,model.rxnECNumbers,...
    model.grRules,model.subSystems,num2cell(model.lb),...
    num2cell(model.ub),num2cell(model.freeEnergy),tags];
xlswrite('Supplementary_Materials.xlsx',A,'Reactions');

% Compile the metabolites information
A = [{'ID','Name','KEGG ID','Formula','Charge'};...
    model.mets,model.metNames,model.metKEGGID,model.metFormulas,...
    num2cell(model.metCharge)];
xlswrite('Supplementary_Materials.xlsx',A,'Metabolites');

% Run an FVA on the H2-consuming model
[minFlux,maxFlux] = fluxVariability(model,100);
% Compile the answers into a sheet with reactions
A = [{'Reaction ID','Minimum Flux','Maximum Flux'};...
    model.rxns,num2cell(minFlux),num2cell(maxFlux)];
xlswrite('Supplementary_Materials.xlsx',A,'FVA on H2');

% Switch to formate and repeat the FVA step
model = switchToFormate(model);
[minFlux,maxFlux] = fluxVariability(model,100);
% Compile the answers into a sheet with reactions
A = [{'Reaction ID','Minimum Flux','Maximum Flux'};...
    model.rxns,num2cell(minFlux),num2cell(maxFlux)];
xlswrite('Supplementary_Materials.xlsx',A,'FVA on Formate');


end

function tags=createRxnConfidenceSheet(model)
%
%Create a list of tags that includes whether each reaction was in the
%original reconstruction, gap-filled to make the original model, or
%manually added later on, plus tags for transport and exchanges
%Load the original model (M_mar)
load('original_model.mat')

%Remove the non-genes from the original model
M_mar=removeGene(M_mar,'Unknown');
M_mar=removeGene(M_mar,'fig');

%Find the reactions that had genes in the orignal model
recon_rxns = setdiff(M_mar.rxns,findRxnsWOGenes(M_mar));

%Find exchange reactions, give them their own tags in both
%models
exc_rxns = model.rxns(findExcRxns(model));

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
% Use "Kbase" to indicate it's a gene from Kbase
[~,idx] = intersect(model.rxns,recon_rxns);
for i=1:length(idx)
    tags{idx(i)}='KBase';
end

%Use 'EX' for 'Exchange Reaction'
% Use "Exchange"
[~,idx] = intersect(model.rxns,exc_rxns);
for i=1:length(idx)
    tags{idx(i)}='Exchange';
end

% Use 'GF' for 'Gapfilled Reaction'
% Say "Gapfill"
[~,idx] = intersect(model.rxns,gapfill_rxns);
for i=1:length(idx)
    tags{idx(i)}='Gapfill';
end

%Step 3: Basically add "AM" for reactions that were "Added Manually"
% Just write "Manual Addition"
[~,idx] = setdiff(model.rxns,M_mar.rxns);
for i=1:length(idx)
    tags{idx(i)}='Manual Addition';
end

%Also do it for anything that's in the manual set
[~,idx] = intersect(model.rxns,manual_rxns);
for i=1:length(idx)
    tags{idx(i)}='Manual Addition';
end


%%%
% I want to do something different for transport: if they're Kbase things,
% then leave that tag; if they're in the original as Gapfills, call them
% "physiological"; if they're manually added, call them "Manual Addition"
% Label transport reactions as Transport for now
trans_rxns = findTransRxnsMOD(model);
[~,idx] = intersect(model.rxns,intersect(trans_rxns,gapfill_rxns));
for i=1:length(idx)
    tags{idx(i)}='Physiological';
end

end