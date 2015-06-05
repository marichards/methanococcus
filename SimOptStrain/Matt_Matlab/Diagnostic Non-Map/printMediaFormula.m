function [media_rxns,amounts] = printMediaFormula(model,printFlag)

%Print out the media formulation the model is currently growing on when
%given the model
%
%Inputs
%model: a COBRA model structure
%
%Outputs
%media_rxns: a cell array containing the names of all reactions that can be
%considered in the "media"
%amounts: an array with reaction lower bounds corresponding to the
%reactions in media_rxns

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

if printFlag
    fprintf('\nPossible Media Reactions\n\n')
    for i=1:length(media_rxns)
        fprintf('%s\t%f\n',media_rxns{i},amounts(i))
    end
    fprintf('\n')
end
    