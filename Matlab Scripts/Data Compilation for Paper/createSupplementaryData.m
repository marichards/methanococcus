function createSupplementaryData(model)
%%
% Takes in the M. maripaludis model and creates the automated portion of
% the Supplementary Materials Excel workbook containing: 
% Reactions: info for all reactions, including cross references, formulas,
% GPRs, subsystems, names, bounds,free energy (if applicable)
% Metabolites: info for all metabolites, including cross references,
% formulas, charges
% FVA Results on H2: bounds for reaction fluxes on H2+CO2 media
% FVA Results on Formate: bounds for reaction fluxes on formate media
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% Matthew Richards, 09/28/2015


% Compile the reactions information
formulas = printRxnFormula(model,model.rxns,false,false,true);
% Add the origin confidence data from other code (below)
tags = createRxnConfidenceSheet(model);
% Add the reaction probabilities from other code (below)
probs = findRxnProbabilities(model);
A = [{'ID','Name','KEGG ID','Formula','EC Number(s)','GPR Rules',...
    'Subsystem','Lower Bound','Upper Bound','Free Energy (kJ/mmol)',...
    'Likelihood Score','Origin Tag'};...
    model.rxns,model.rxnNames,model.rxnKEGGID,formulas,model.rxnECNumbers,...
    model.grRules,model.subSystems,num2cell(model.lb),...
    num2cell(model.ub),num2cell(model.freeEnergy),probs,tags];
xlswrite('Supplementary_Materials.xlsx',A,'Reactions');

% Compile the metabolites information
A = [{'ID','Name','KEGG ID','ChEBI ID','Formula','Charge'};...
    model.mets,model.metNames,model.metKEGGID,model.metChEBIID,...
    model.metFormulas,num2cell(model.metCharge)];
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

function probs = findRxnProbabilities(model)
%%
% Create a list of probabilities that specify the likelihood of each 
% reaction in the rxn probs object from Kbase. Print N/A for anything
% without complexes
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% OUTPUT
% probs: a list of rxn probability tags for each reaction in the 
% M. maripaludis model that indicates its probability (if it exists) and
% specifies "N/A" otherwise
%
%
% Matthew Richards, 09/28/2015

% Load the probabilities data
load('all_rxn_probs.mat')

% Create a list of probabilities
probs = cell(size(model.rxns));

% Loop through reactions
for i = 1:length(model.rxns)
    % Look for the reaction in the rxn probs object
    if ismember(model.rxns{i},all_reactions)
        % If it's there, then grab its probability and convert it to a
        % string
        [~,idx] = intersect(all_reactions,model.rxns{i});
        probs{i} = num2str(probability1(idx));
    else 
        % If not, it's N/A because we can't find the complexes
        probs{i} = 'N/A';
    end
end
end

function tags=createRxnConfidenceSheet(model)
%%
% Create a list of tags that includes whether each reaction was in the
% original reconstruction, gap-filled to make the original model, or
% manually added later on, plus tags for transport and exchanges
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% OUTPUT
% tags: a list of origin tags for each reaction in the M. maripaludis model
% that indicates the origin or every reaction
%
% Matthew Richards, 09/28/2015


%Load the original draft model (M_mar)
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