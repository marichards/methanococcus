function formatRxnsForCytoscape(model,rxns,filename)

%Format a set of reactions for cytoscape
%Want a text file with clearly-delimited reactant/reaction or
%reaction/product

%Input:
%rxns - a cell array containing reaction IDs for a COBRA model
%filename - a character string used to name the output file

%Output: 
%A text file named 'filename.txt' that contains the selected reactions
%formatted for Cytoscape visualization

%Open the text file
file_id = fopen(sprintf('%s.txt',filename),'w');

%First loop through the reactions
for i =1:length(rxns)
   
    %Find all the metabolites that are reactants
    reactants = findRxntsFromRxn(model,rxns{i});
    
    for j=1:length(reactants)
    %Print all those reactants
    fprintf(file_id,'%s\t%s\n',reactants{j},rxns{i});
    end
    
    %Find all the metabolites that are products
    products = findProdsFromRxn(model,rxns{i});
    
    for j=1:length(products)
        %Print all those products
        fprintf(file_id,'%s\t%s\n',rxns{i},products{j});
    end
end

end

function reactants = findRxntsFromRxn(model,rxn)

%Find the metabolites in a given reaction that are reactants

%First find the reaction index
[~,idx] = intersect(model.rxns,rxn);

%Using that index, find all indices of reactants (less than 1)
reactants = model.mets(model.S(:,idx)<0);

end

function products = findProdsFromRxn(model,rxn)

%Find the metabolites in a given reaction that are reactants

%First find the reaction index
[~,idx] = intersect(model.rxns,rxn);

%Using that index, find all indices of products (greater than 1)
products = model.mets(model.S(:,idx)>0);

end

