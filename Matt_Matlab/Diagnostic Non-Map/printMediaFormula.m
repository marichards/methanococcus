function [media_rxns,amounts] = printMediaFormula(model,printFlag)

% Print out the media formulation the model is currently growing on when
% given the model. This will return the possible media components, which
% are those with lower bounds less than zero. 
%
% INPUT
% model: a COBRA Toolbox model structure
%
% OPTIONAL INPUT:
% printFlag: a binary variable that, when set to true, will cause the
% function to print out the available media components. (Default = true)
%
% OUTPUT
% media_rxns: a cell array containing the names of all reactions that can
% be considered in the "media"
% amounts: an array with reaction lower bounds corresponding to the
% reactions in media_rxns
%
% Matthew Richards, 09/24/2015

% Check for the print flag, set it to true if there's nothing specified
if (nargin <2)
    printFlag = 1;
end

%Find all the exchange reactions
exc_rxns = model.rxnNames(findExcRxns(model));

%Find all the indices
[exc_rxns,idx] = intersect(model.rxnNames,exc_rxns);

%Initiate media_rxns
media_rxns = {};
amounts = [];
%Run through and save things that can have uptake
for i = 1:length(exc_rxns)
    %Check if LB is less than 0
    if model.lb(idx(i))<0
        %Then save it and the amount
        media_rxns = [media_rxns;exc_rxns{i}];
        amounts = [amounts;model.lb(idx(i))];
  
    end
end

% If instructed, print out the media reactions
if printFlag
    fprintf('\nPossible Media Reactions\n\n')
    for i=1:length(media_rxns)
        fprintf('%s\t%f\n',media_rxns{i},amounts(i))
    end
    fprintf('\n')
end
    