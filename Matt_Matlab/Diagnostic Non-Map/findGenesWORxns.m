function genes = findGenesWORxns(model)

genes = {};

for i =1:length(model.genes)
   
    %If there are no genes
    if sum(model.rxnGeneMat(:,i)) == 0
        
        %Then save them
        genes = [genes;model.genes{i}]; 
    end
    

end
