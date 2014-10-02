function rxns = findRxnsWOGenes(model)

rxns = {};

for i =1:length(model.rxns)
   
    %If there are no genes
    if isempty(model.grRules{i})
        
        %Then save them
        rxns = [rxns;model.rxns{i}]; 
    end
    

end

%Pull out exchanges
rxns = setdiff(rxns,model.rxns(findExcRxns(model)));