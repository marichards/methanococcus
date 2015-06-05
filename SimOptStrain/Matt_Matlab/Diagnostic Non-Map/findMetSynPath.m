function [all_rxns,rxnNames,fluxes] = findMetSynPath(model,met,solution,excluded,threshold)

%Select any metabolite in the model and pull the highest-flux pathways for
%that metabolite, then use those to trace where the metabolite is coming
%from and going to

%Inputs
%model: a COBRA model structure
%met: the metabolite of interest
%solution: flux solution from running optimizeCbModel

%Outputs
%all_rxns: List of reactions making/using a metabolite
%fluxes: List of fluxes through those reactions

%Find the exchange reaction IDs. These will be the end points of a
%"pathway"
%May not need these
%exchanges = model.rxns(findExcRxns(model,'True'));

%If no excluded reactions, excluded is empty
if (nargin < 4)
    excluded={''};
end

if (nargin < 5)
    threshold = 1;
end



%List the cofactor metabolites so we can avoid them
%Pull out things with frequencies higher than 30 reactions
cofactors = {};
for i=1:length(model.mets)
freq=length(findRxnsFromMets(model,model.mets{i}));
if freq>20
cofactors=[cofactors;model.mets{i}];
end
end

%If the metabolite is in the cofactors, return an error
if ismember(met,cofactors)
    
    error('Please supply a primary metabolite');
end

%Start with the metabolite...make it a cell array
mets = {met};

%Create all_rxns array and all_mets array
all_rxns={};
all_mets={met};
%Now enter a while loop...for now, make it while its smaller than the whole
%model, so it'll never be violated
while length(all_rxns)<20
    
    %%Upgrade using metMassBalance code
    %Use metabolites to find reactions with notable fluxes using metMassBalance
    %Initiate the array
    %%%Change on 1/15/2015: Only find things that make it(?)
    rxns = {};
    for i=1:length(mets)
        %Narrow to only reactions that PRODUCE it
        [new_rxns,amounts] = metMassBalance(model,mets(i),solution,threshold);
        rxns = [rxns; new_rxns(amounts>0)];
    end
    
    
    
    %Narrow it down to only new reactions
    rxns = setdiff(rxns,all_rxns);
    %Pull out excluded reactions
    rxns = setdiff(rxns,excluded);
    
    %If there are no new ones, break the loop
    if length(rxns)<1
        break
    end

    %Add to the all_rxns
    all_rxns = [all_rxns;rxns];
  
    %Use the new reactions to find metabolites
    mets = findMetsFromRxns(model,rxns);
    
    %Make sure to pull out the cofactor metabolites
    mets = setdiff(mets,cofactors);
    
    %Narrow it down to only new metabolites
    mets = setdiff(mets,all_mets);
    
    %If there are no new ones, break the loop
    if length(mets)<1
        break
    end
    
    %Add to the all_mets
    all_mets = [all_mets;mets];
        
end

%Pull out the fluxes
[all_rxns,all_idx] = intersect(model.rxns,all_rxns);
fluxes=solution.x(all_idx);
rxnNames = model.rxnNames(all_idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


