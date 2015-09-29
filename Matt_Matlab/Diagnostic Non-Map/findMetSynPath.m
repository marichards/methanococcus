function [all_rxns,rxnNames,fluxes] = findMetSynPath(model,met,solution,excluded,threshold,max_rxns)

% Select any metabolite in the model and pull the highest-flux pathways for
% that metabolite, then use those to trace where the metabolite is coming
% from and going to in the model
%
% INPUT:
% model: a COBRA Toolbox model structure
% met: a metabolite of interest in the supplied model
% solution: a flux distribution solution for the supplied model
%
% OPTIONAL INPUT: 
% excluded: a list of reactions to excluded from the metabolite synthesis
% pathway. Often, this includes the biomass reaction to avoid including the
% many metabolites outlined in that reaction. (Default = {''}) 
% threshold: a value scaled from 0-1 that serves as the cutoff for
% following a flux. This value sets the minimum percentage of flux to
% follow in synthesis. For example, setting the value to "0.5" will follow
% all non-excluded reactions that carry at least 50% of the flux of a
% primary metabolite from reactions already in the pathway, but will not
% follow any reactions below that threshold. (Default = 0.1)
% max_rxns: a maximum  value on the number of reactions that can be
% included in the pathway. The function will discontinue tracing the
% synthesis path as soon as this number is exceeded (Default = 20)
%
% OUTPUT: 
% all_rxns: list of reactions making/using a metabolite
% rxnNames: list of names corresponding to the reactions in all_rxns
% fluxes: list of fluxes through the reactions in all_rxns
%
% Matthew Richards, 09/24/2015


% If no excluded reactions, excluded is empty
if (nargin < 4)
    excluded={''};
end
% If no threshold, set it at 10%
if (nargin < 5)
    threshold = 0.1;
end
% If no maximum reactions given, set it at 20
if (nargin <6)
    max_rxns = 20;
end

% List the cofactor metabolites so we can avoid them
% Pull out things with frequencies higher than 30 reactions
cofactors = {};
for i=1:length(model.mets)
freq=length(findRxnsFromMets(model,model.mets{i}));
if freq>20
cofactors=[cofactors;model.mets{i}];
end
end

% If the metabolite is in the cofactors, return an error
if ismember(met,cofactors)
    
    error('Please supply a primary metabolite');
end

% Start with the metabolite...make it a cell array
mets = {met};

% Create all_rxns array and all_mets array
all_rxns={};
all_mets={met};
% Now enter a while loop...for now, make it while its smaller than the whole
% model, so it'll never be violated
while length(all_rxns)<max_rxns
    
    % Use metabolites to find reactions with notable fluxes using metMassBalance
    % Initiate the array
    %%% Change on 1/15/2015: Only find pathways that synthesize it
    rxns = {};
    for i=1:length(mets)
        % Narrow to only reactions that PRODUCE it
        [new_rxns,amounts] = metMassBalance(model,mets(i),solution,threshold);
        rxns = [rxns; new_rxns(amounts>0)];
    end
    
    % Narrow it down to only new reactions
    rxns = setdiff(rxns,all_rxns);
    % Pull out excluded reactions
    rxns = setdiff(rxns,excluded);
    
    % If there are no new ones, break the loop
    if length(rxns)<1
        break
    end

    % Add to the all_rxns
    all_rxns = [all_rxns;rxns];
  
    % Use the new reactions to find metabolites
    mets = findMetsFromRxns(model,rxns);
    
    % Make sure to pull out the cofactor metabolites
    mets = setdiff(mets,cofactors);
    
    % Narrow it down to only new metabolites
    mets = setdiff(mets,all_mets);
    
    % If there are no new ones, break the loop
    if length(mets)<1
        break
    end
    
    % Add to the all_mets
    all_mets = [all_mets;mets];
        
end

% Pull out the fluxes
[all_rxns,all_idx] = intersect(model.rxns,all_rxns);
fluxes=solution.x(all_idx);
rxnNames = model.rxnNames(all_idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


