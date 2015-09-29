function [essentialRxns,essentialDirections] =findEssentialRxnsFromMet(model,met,excludedRxns)

% Takes in a metabolite and a model, finds all reactions for that
% metabolite, then uses the diagnoseMissingMets code to narrow it to all
% reactions that are essential for that metabolite's production or
% consumption.  Returns the list of those reactions and whether or not
% they're a sink or a source
%
% INPUT
% model: a COBRA Toolbox model structure
% met: a metabolite of interest
%
% OPTIONAL INPUT
% excludedRxns: a list of reactions to exclude from the list of essential
% reactions. The intent of this input is to exclude model artifacts, such
% as the biomass equation, and focus on essential internal reactions
% (Default = {})
%
% OUTPUT
% essentialRxns: a list of reactions essential for production or
% consumption of the metabolite of interest
% essentialDirections: designation of 'source' or 'sink' corresponding to 
% each reaction
%
% Matthew Richards, 09/24/2015

% Step 1: Find all reactions for met and pull out excluded reactions if
% they exist
rxns = findRxnsFromMets(model,met);

if (nargin)>2
    rxns = setdiff(rxns,excludedRxns);
end

% Step 2: Create arrays to collect essential reactions and directions
essentialRxns={};
essentialDirections = {};

% Step 3: Cycle through these reactions
for i = 1:length(rxns)
    
    % Step 4: Check for the missing metabolites in each
    [missingMets,directions] = diagnoseMissingMet(model,rxns{i});
    
    % Step 5: If the met is in the missingMets...
    if ismember(met,missingMets)
        
       % Step 6: Pull out the index
       [~,idx] = intersect(missingMets,met);
       
       % Step 7: Store the reaction and direction in the arrays
       essentialRxns=[essentialRxns; rxns{i}];
       essentialDirections = [essentialDirections; directions{idx}];            
    end    
end