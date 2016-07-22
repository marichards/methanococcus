function [transRxns,nonTransRxns] = findTransRxnsMOD(model,inclExc)
% findTransRxnsMOD identifies all transport reactions in the M. maripaludis
% model, or any other model that includes compartments 'c0','p0','e0'. 
%
% [transRxns,nonTransRxns] = findTransRxns(model,inclExc)
%
%INPUT
% model             COBRA model structure
%
%OPTIONAL INPUT
% inclExc           include exchange reactions as transport?
%                   (Default = false)
%
%OUTPUT
% transRxns         all transport reactions in the model
% nonTransRxns      all non-transport reactions in the model
%
% right now, this function only works with models the compartments [c0],
% [p0], and [e0]
%
% Jeff Orth  8/31/07
%
% Modified by Matthew Richards to fit compartments for M. maripaludis,
% 09/28/2015

if nargin < 2
    inclExc = false;
end

isTrans = zeros(1,length(model.rxns));

for i = 1:length(model.rxns)
    mets = model.mets(find(model.S(:,i)));
    cMets = regexp(mets,'\[c0\]');
    hasCs = ~isempty([cMets{:}]);
    pMets = regexp(mets,'\[p0\]');
    hasPs = ~isempty([pMets{:}]);
    eMets = regexp(mets,'\[e0\]');
    hasEs = ~isempty([eMets{:}]);
    
    if (sum([hasCs,hasPs,hasEs]) > 1) || hasEs&&inclExc
        isTrans(i) = 1;
    end
end

transRxns = model.rxns(isTrans==1);
nonTransRxns = model.rxns(isTrans==0);



