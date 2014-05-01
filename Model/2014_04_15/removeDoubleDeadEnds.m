function [new_model,removedRxns,removedMets]=removeDoubleDeadEnds(model)

%%%%%%%%
%Find dead end metabolites
deadMets=model.mets(detectDeadEnds(model));

%Make arrays for reactions
pos_rxns={};
neg_rxns={};

%For each metabolite
for i=1:length(deadMets)
%Find the proper index in the model
[~,idx]=intersect(model.mets,deadMets{i});

%Find reaction that it's in (there will be only one)
rxn=find(model.S(idx,:));
%Add to the proper array by testing if it's positive or negative in the S
%matrix
if model.S(idx,rxn)<0
    neg_rxns = [neg_rxns; model.rxns{rxn}];
else
    
    pos_rxns = [pos_rxns; model.rxns{rxn}];
end
end

%Find the reactions that are in both and remove them
removedRxns = intersect(pos_rxns,neg_rxns);

%Remove the reactions!
new_model=removeRxns(model,removedRxns);

%Find the removed metabolites
removedMets=setdiff(model.mets,new_model.mets);

