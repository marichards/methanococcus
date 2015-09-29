function genes = findGenesWORxns(model)

% Searches through a model and finds all genes that lack reactions, or
% "orphan genes". These genes are unconnected to the stoichiometric matrix
% and may be candidates for removal. 
%
% INPUT
% model: a COBRA Toolbox model structure
%
% OUTPUT
% genes: a list of orphan genes in the model that lack reactions
%
% Matthew Richards, 09/24/2015

genes = {};

for i =1:length(model.genes)
   
    %If there are no genes
    if sum(model.rxnGeneMat(:,i)) == 0
        
        %Then save them
        genes = [genes;model.genes{i}]; 
    end
    

end
