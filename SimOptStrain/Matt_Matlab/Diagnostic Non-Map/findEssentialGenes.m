function essentials = findEssentialGenes(model,genes,cutoff)

%Based on a cutoff for percentage of growth considered "lethal", find the
%essential genes in the given list of genes

%Qualify if there isn't a cutoff; make it 0
if (nargin<3)
    cutoff = 0;
end

%If there's no genes given, then just use all of them
if (nargin<2)
    genes = model.genes;
end

% Run the single gene deletion
grRatio = singleGeneDeletion(model,'FBA',genes);

%Pull out things where the growth ratio is below the cutoff
essentials = genes(grRatio<=cutoff);